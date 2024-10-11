use crate::quantized_mesh_tile::QuantizedMeshTile;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

#[derive(PartialEq)]
pub enum CoordType {
    Geocentric,
    Geodetic,
}

pub fn write_csv(qmt: &QuantizedMeshTile, file_path: &Path, coord_type: CoordType) -> io::Result<()> {
    let mut vertices = String::new();
    let mut triangles = String::new();
    let mut vfile = File::create(file_path.with_extension("vertices.csv"))?;
    let mut tfile = File::create(file_path.with_extension("triangles.csv"))?;

    triangles.push_str("v0, v1, v2\n");
    for t in qmt.quantized_mesh.vertex_data.triangle_index.iter() {
        triangles.push_str(&format!("{}, {}, {}\n", t[0], t[1], t[2]));
    }
    
    if coord_type == CoordType::Geocentric {
        vertices.push_str("x, y, z\n");
        for vlatlon in qmt.vertices_as_geodetic_point3() {
            let v = vlatlon.to_ecef(&qmt.ellipsoid);
            vertices.push_str(&format!("{}, {}, {}\n", v.0.x, v.0.y, v.0.z));
        }
    }
    else {
        vertices.push_str("lat, lon, height\n");
        for vertex in qmt.vertices_as_geodetic_point3() {
            vertices.push_str(&format!(
                "{:.12}, {:.12}, {:.12}\n",
                vertex.0.x, vertex.0.y, vertex.0.z
            ));
        }
    }
    
    let _ = vfile.write_all(vertices.as_bytes());
    let _ = tfile.write_all(triangles.as_bytes());
    Ok(())
}
