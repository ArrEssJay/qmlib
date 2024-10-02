use std::{error::Error, fs::File, path::Path};
use tiff::{
    encoder::{colortype::Gray32Float, TiffEncoder},
    tags::Tag,
};

use crate::{
    geometry::GeodeticLatLon, raster::{raster_size_pixels, rasterise}, quantized_mesh_tile::{QuantizedMeshTile, CRS}
};

pub fn write_tiff(
    qmt: &QuantizedMeshTile,
    filename: &Path,
    scale_shift: u16,
) -> Result<(), Box<dyn Error>> {
    let raster_size = raster_size_pixels(scale_shift);
    let raster = rasterise(qmt, scale_shift);

    let mut tiff: TiffEncoder<File> = TiffEncoder::new(File::create(filename)?)?;
    let mut image = tiff
        .new_image::<Gray32Float>(raster_size.into(), raster_size.into())
        .unwrap();

    // TODO - Handle Projected Web Mercator
    let geo_key_directory = if qmt.crs == CRS::Epsg4326 {
     vec![
    1, 1, 0, 2,          // Header: Version 1.1, 8 keys
    1024, 0, 1, 2,      // GTRasterTypeGeoKey 2 - geographic
    2048, 0, 1, 4979,   // GeodeticCRSGeoKey, value 4979 (WGS 84 3D)
    4096, 0, 1, 4979,   // VerticalGeoKey, value 4979 (WGS 84 3D)
    ]} else {
        panic!("Unsupported CRS: {:?}", qmt.crs);
    };

    // Now define the bounding rectangle values
    let bounding_rectangle = &qmt.bounding_rectangle; // Assuming this is your calculated bounding box
    let ll = bounding_rectangle.lower_left.to_degrees();
    let ur = bounding_rectangle.upper_right.to_degrees();

    let min_x = ll.lat();
    let min_y = ll.lon();
    let max_x = ur.lat();
    let max_y = ur.lon();
 
    // 0,0 is UL
    // crs axis order is lat,lon,alt
    // Define ModelTiepoint
    // Technically we only need one of these
    let tiepoints: Vec<f64> = vec![
        0.0,    // UL UV x
        0.0,    // UL UV y
        0.0,    // UL UV Z
        *min_x, // UL world x
        *max_y, // UL world y
        0.0,    // UL world z
        (raster_size as f64), //LR UV x
        (raster_size as f64), //LR UV y
        0.0,  //LR UV z
        *max_x,  // LR world x
        *min_y,  // LR world y
        0.0,     // LR world z
    ];

    // Define ModelPixelScale (assuming square pixels)
    let pixel_size_x = (*max_x - *min_x) / raster_size as f64; // Calculate pixel size in X
    let pixel_size_y = (*max_y - *min_y) / raster_size as f64; // Calculate pixel size in Y
    let pixel_scale = vec![pixel_size_x, pixel_size_y, 1.0]; // Set Z pixel size to 1 as Z values require no scaling

    // Debugging: Print bounding rectangle coordinates
    println!("Bounding Rectangle:");
    println!("Lower Left: ({}, {})", min_x, min_y);
    println!("Upper Right: ({}, {})", max_x, max_y);
    println!("Scale: ({:?})", pixel_scale);
    println!("raster_size: ({:?})", raster_size);

    image.encoder().write_tag(Tag::ModelPixelScaleTag, &pixel_scale[..]).unwrap();

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
    let expected_size = (raster_size as u32) * (raster_size as u32);
    if raster.len() as u32 != expected_size {
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
