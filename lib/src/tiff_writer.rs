use std::{error::Error, fs::File, path::Path};
use tiff::{
    encoder::{colortype::Gray32Float, TiffEncoder},
    tags::Tag,
};

use crate::quantized_mesh_tile::{QuantizedMeshTile, CRS};

use compute_shader::RasterParameters;
use compute_shader_interface::VertexBuffers;
use compute_shader_interface::rasterise;
pub use compute_shader_interface::Rasteriser; // re-export Rasteriser for a unified interface

pub fn write_tiff(
    qmt: &QuantizedMeshTile,
    filename: &Path,
    raster_scale_factor: u32,
    rasteriser: Rasteriser
) -> Result<(), Box<dyn Error>> {
    let u = &qmt
        .quantized_mesh
        .vertex_data
        .u
        .iter()
        .map(|&x| x as u32)
        .collect::<Vec<u32>>();
    let v = &qmt
        .quantized_mesh
        .vertex_data
        .v
        .iter()
        .map(|&x| x as u32)
        .collect::<Vec<u32>>();
    let attribute = &qmt
        .quantized_mesh
        .vertex_data
        .height
        .iter()
        .map(|&x| x as u32)
        .collect::<Vec<u32>>();
    let indices = &qmt
        .quantized_mesh
        .vertex_data
        .triangle_index
        .iter()
        .flatten().copied()
        .collect::<Vec<u32>>();



    let vertex_buffers =  VertexBuffers {
         u,
         v,
         attribute,
         indices,
    };

    let params = RasterParameters::new(
        raster_scale_factor,
        32767,
        qmt.quantized_mesh.header.min_height,
        qmt.quantized_mesh.header.max_height,
        32767, // max height in quantised mesh
        qmt.quantized_mesh.vertex_data.vertex_count,
        qmt.quantized_mesh.vertex_data.triangle_count,
    );
    let raster_size = params.scaled_raster_size();


    println!("Raster Parameters: {:?}", params);
    println!("Raster Scaled Size: {:?}x{:?}", raster_size, raster_size);

    let raster = rasterise(vertex_buffers, &params, rasteriser);
    let mut tiff: TiffEncoder<File> = TiffEncoder::new(File::create(filename)?)?;
    let mut image = tiff
        .new_image::<Gray32Float>(raster_size, raster_size)
        .unwrap();

    // TODO - Handle Projected Web Mercator
    let geo_key_directory = if qmt.crs == CRS::Epsg4326 {
        vec![
            1, 1, 2, 5, // Header: Version 1.2, 4 keys
            1024, 0, 1, 2, // GTModelTypeGeoKey = 2 (geographic unprojected 2D)
            2048, 0, 1, 4326, // GeodeticCRSGeoKey = 4326 (WGS 84 2D)
            4096, 0, 1,
            4979, // VerticalGeoKey = 4979 (WGS 84 3D) (used here to document use of ellipsoidal height)
            4099, 0, 1, 9001, // VerticalUnitsGeoKey = 9001 (Linear_Meter)
            1025, 0, 1, 2, // GTRasterTypeGeoKey = 2 (RasterPixelIsPoint)
        ]
    } else {
        panic!("Unsupported CRS: {:?}", qmt.crs);
    };

    // Mapping raster space to world space. From the spec:
    // u=>latitude.  0=>lat_min (West), 32767=>lat_max (East)
    // v=>longitude. 0=>lon_min (South), 32767=>lon_max (North)

    // https://docs.ogc.org/is/19-008r4/19-008r4.html#_requirements_class_modelpixelscaletag
    //
    // A positive ScaleX in the ModelPixelScaleTag SHALL indicate that model space X coordinates
    // increase as raster space I indices increase. This is the standard horizontal relationship
    // between raster space and model space. A positive ScaleY in the ModelPixelScaleTag SHALL
    // indicate that model space Y coordinates decrease as raster space J indices increase.
    // This is the standard vertical relationship between raster space and model space.
    // The ScaleZ is primarily used to map the pixel value of a digital elevation model
    // into the correct Z-scale (in other words a Z-Scaling factor).

    // Therefore we require:
    // y scale: negative
    // x scale: positive
    // z: 1 (float32 height 1:1 maps to z)

    // PixelIsPoint
    // u:  When the u value is 0, the vertex is on the Western edge of the tile.
    // If a point-pixel image were to be displayed on a display device with pixel cells having
    // the same size as the raster spacing, then the upper-left corner of the displayed image
    // would be located in raster space at (-0.5, -0.5).

    let bounding_rectangle = &qmt.bounding_rectangle;
    let ll = bounding_rectangle.lower_left.to_degrees();
    let ur = bounding_rectangle.upper_right.to_degrees();

    let min_lat = &ll.0;
    let min_lon = &ll.1;
    let max_lat = &ur.0;
    let max_lon = &ur.1;

    // ModelTiepoint
    let tiepoints: Vec<f64> = vec![
        0.0,      //LL UV x
        0.0,      //LL UV y
        0.0,      //LL UV z
        *min_lon, // LL world x
        *max_lat, // LL world y
        0.0,      // LL world z
    ];

    // Define ModelPixelScale (assuming square pixels)
    let height = max_lat - min_lat;
    let width = max_lon - min_lon;
    let pixel_size_x = width / f64::from(raster_size); // Calculate pixel size in X
    let pixel_size_y = height / f64::from(raster_size); // Calculate pixel size in Y
    let pixel_scale = [pixel_size_x, pixel_size_y, 1.0]; // Set Z pixel size to 1 as Z values require no scaling

    #[cfg(debug_assertions)]
    {
        println!("Tie Points: {tiepoints:?}");
        println!("Bounding Rectange: {ll:?} {ur:?}");
        println!("Bounding Rectangle Height: {height:?}");
        println!("Bounding Rectangle Width: {width:?}");

        println!("Scale: ({pixel_scale:?})");
        println!("raster_size: ({raster_size:?})");
    }

    image
        .encoder()
        .write_tag(Tag::ModelPixelScaleTag, &pixel_scale[..])
        .unwrap();

    // Write the ModelTiepointTag
    image
        .encoder()
        .write_tag(Tag::ModelTiepointTag, &tiepoints[..])
        .unwrap();

    image
        .encoder()
        .write_tag(Tag::GeoKeyDirectoryTag, &geo_key_directory[..])
        .unwrap();

    // Ensure we have the correct dimensions for the raster
    let expected_size = raster_size * raster_size;
    if u32::try_from(raster.len()) != Ok(expected_size) {
        return Err(format!(
            "Raster size does not match expected dimensions: expected {}, found {}",
            expected_size,
            raster.len()
        )
        .into());
    }
    //image.rows_per_strip(2).unwrap();

    // Write the image data as Gray32Float
    image.write_data(&raster)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use compute_shader_interface::Rasteriser;

    use crate::quantized_mesh_tile::load_quantized_mesh_tile;
    use crate::{test_utils::test_data::qmt_test_chess, tiff_writer::write_tiff};
    use crate::twodm_writer::write_2dm;
    #[test]
    fn test_write_tiff() {
        let mesh = qmt_test_chess();

        let raster_scale_factor: u32 = 5; // Example raster_scale_factor for rasterisation
        let result = write_tiff(&mesh, &PathBuf::from("./test_cpu.tiff"), raster_scale_factor, Rasteriser::CPU);
        assert!(result.is_ok(), "CPU TIFF generation failed: {result:?}");

        let result = write_tiff(&mesh, &PathBuf::from("./test_gpu.tiff"), raster_scale_factor, Rasteriser::GPU);
        assert!(result.is_ok(), "GPU TIFF generation failed: {result:?}");

        let result = write_2dm(&mesh, &PathBuf::from("./test.2dm"));
        assert!(result.is_ok(), "2dm generation failed: {result:?}");
    } 

    #[test]

    fn test_write_tiff_terrain() {
        let path = PathBuf::from(format!(
            "{}/../test/terrain_data/a/15/59489/9692.terrain",
            env!("CARGO_MANIFEST_DIR")
        ));
        let mesh = load_quantized_mesh_tile(&path).unwrap();

        let raster_scale_factor: u32 = 5; // Example raster_scale_factor for rasterisation
        let result = write_tiff(&mesh, &path, raster_scale_factor, Rasteriser::CPU);

        assert!(result.is_ok(), "TIFF generation failed: {result:?}");
    }
}
