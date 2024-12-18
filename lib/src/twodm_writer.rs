use std::fs::File;
use std::io::{self, Write};
use std::path::Path;
use crate::quantized_mesh_tile::{QuantizedMeshTile, CRS};

pub fn write_2dm(tile: &QuantizedMeshTile, file_path: &Path) -> io::Result<()> {
    // Write the .2dm file
    let mut file = File::create(file_path)?;
    
    // Write header
    writeln!(file, "MESH2D")?;
    
    // Use the existing vertices in degrees
    let vertices = tile.vertices_as_geodetic_point3();
    
    // Write nodes (coordinates in degrees)
    for (i, vertex) in vertices.iter().enumerate() {
        writeln!(
            file,
            "ND {} {} {} {}",
            i + 1,
            vertex.lon().to_degrees(),
            vertex.lat().to_degrees(),
            vertex.height()
        )?;
    }
    
    // Write elements (triangles)
    let index_data = &tile.quantized_mesh.vertex_data.triangle_index;
    for (i, indices) in index_data.iter().enumerate() {
        writeln!(
            file,
            "E3T {} {} {} {}",
            i + 1,
            indices[0] + 1,
            indices[1] + 1,
            indices[2] + 1
        )?;
    }
    
    // Write the .prj file with CRS information
    let prj_path = file_path.with_extension("prj");
    let mut prj_file = File::create(prj_path)?;
    
    // Select the appropriate WKT based on the CRS
    let crs_wkt = match tile.crs {
        CRS::Epsg4326 => {
            // WKT for EPSG:4326 (WGS84)
            r#"GEOGCS["WGS 84",
    DATUM["WGS_1984",
        SPHEROID["WGS 84",6378137,298.257223563,
            AUTHORITY["EPSG","7030"]],
        AUTHORITY["EPSG","6326"]],
    PRIMEM["Greenwich",0,
        AUTHORITY["EPSG","8901"]],
    UNIT["degree",0.0174532925199433,
        AUTHORITY["EPSG","9122"]],
    AUTHORITY["EPSG","4326"]]"#
        },
        CRS::Epsg3857 => {
            // WKT for EPSG:3857 (Web Mercator)
            r#"PROJCS["WGS 84 / Pseudo-Mercator",
    GEOGCS["WGS 84",
        DATUM["WGS_1984",
            SPHEROID["WGS 84",6378137,298.257223563,
                AUTHORITY["EPSG","7030"]],
            AUTHORITY["EPSG","6326"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.0174532925199433,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4326"]],
    PROJECTION["Mercator_1SP"],
    PARAMETER["central_meridian",0],
    PARAMETER["scale_factor",1],
    PARAMETER["false_easting",0],
    PARAMETER["false_northing",0],
    UNIT["metre",1,
        AUTHORITY["EPSG","9001"]],
    AXIS["X",EAST],
    AXIS["Y",NORTH],
    AUTHORITY["EPSG","3857"]]"#
        },
    };
    
    prj_file.write_all(crs_wkt.as_bytes())?;
    
    Ok(())
}