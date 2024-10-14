use std::{error::Error, fs::File, path::Path};
use tiff::{
    encoder::{colortype::Gray32Float, TiffEncoder},
    tags::Tag,
};

use crate::{
    interpolator::InterpolationMethod,
    quantized_mesh_tile::{QuantizedMeshTile, CRS},
    raster::{raster_dim_pixels, rasterise},
};

pub fn write_tiff(
    qmt: &QuantizedMeshTile,
    filename: &Path,
    scale_shift: u16,
    method: InterpolationMethod,
) -> Result<(), Box<dyn Error>> {
    let raster_size = raster_dim_pixels(scale_shift);

    let raster = rasterise(qmt, scale_shift, &method);
    let mut tiff: TiffEncoder<File> = TiffEncoder::new(File::create(filename)?)?;
    let mut image = tiff
        .new_image::<Gray32Float>(raster_size.into(), raster_size.into())
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
        *min_lat, // LL world x
        *max_lon, // LL world y
        0.0,      // LL world z
    ];

    // Define ModelPixelScale (assuming square pixels)
    let height = max_lat - min_lat;
    let width = max_lon - min_lon;
    let pixel_size_x = height / f64::from(raster_size); // Calculate pixel size in X
    let pixel_size_y = width / f64::from(raster_size); // Calculate pixel size in Y
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
    let expected_size = u32::from(raster_size) * u32::from(raster_size);
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

    use crate::{
        interpolator::{self},
        test_utils::test_data::qmt_test_chess,
        tiff_writer::write_tiff,
    };

    #[test]
    fn test_write_tiff() {
        let mesh = qmt_test_chess();
        let path: PathBuf = PathBuf::from("./test.tiff");

        let scale_shift: u16 = 8; // Example scale_shift for rasterisation
        let result = write_tiff(
            &mesh,
            &path,
            scale_shift,
            interpolator::InterpolationMethod::Edge,
        );

        assert!(result.is_ok(), "TIFF generation failed: {result:?}");
    }
}
