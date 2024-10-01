use nalgebra::Point2;
use std::{error::Error, fs::File, path::Path};
use tiff::{
    encoder::{colortype::Gray32Float, TiffEncoder},
    tags::Tag,
};

use crate::{
    geometry::Triangle, raster::interpolate_height_in_triangle, QuantizedMesh, UV_MAX_U16,
};

pub fn get_raster_size(scale_shift: u16) -> u16 {
    (UV_MAX_U16 + 1) >> scale_shift
}

pub fn rasterize(qm: &QuantizedMesh, scale_shift: u16) -> Vec<Option<f32>> {
    let raster_size = get_raster_size(scale_shift);
    let mut raster: Vec<Option<f32>> = vec![None; raster_size as usize * raster_size as usize];
    let v = &qm.vertex_data;
    let heights = &qm.interpolated_height_vertices();
    for triangle_index in 0..v.triangle_index.len() {
        let triangle_indices = &v.triangle_index[triangle_index];

        let triangle = Triangle {
            vertices: [
                Point2::new(
                    (qm.vertex_data.u[triangle_indices[0] as usize] >> scale_shift) as f64,
                    (qm.vertex_data.v[triangle_indices[0] as usize] >> scale_shift) as f64,
                ),
                Point2::new(
                    (qm.vertex_data.u[triangle_indices[1] as usize] >> scale_shift) as f64,
                    (qm.vertex_data.v[triangle_indices[1] as usize] >> scale_shift) as f64,
                ),
                Point2::new(
                    (qm.vertex_data.u[triangle_indices[2] as usize] >> scale_shift) as f64,
                    (qm.vertex_data.v[triangle_indices[2] as usize] >> scale_shift) as f64,
                ),
            ],
        };

        let tri_heights = [
            heights[triangle_indices[0]],
            heights[triangle_indices[1]],
            heights[triangle_indices[2]],
        ];

        let triangle_bounds = triangle.bounding_rect();

        for px in (triangle_bounds.lower_left.x as u16).min(raster_size - 1)
            ..=(triangle_bounds.upper_right.x as u16).min(raster_size - 1)
        {
            for py in (triangle_bounds.lower_left.y as u16).min(raster_size - 1)
                ..=(triangle_bounds.upper_right.y as u16).min(raster_size - 1)
            {
                let point = Point2::new(px as f64, py as f64);
                let raster_idx = (py as usize * raster_size as usize) + px as usize;
                if raster[raster_idx] == None {
                    if let Some(interpolated_height) =
                        interpolate_height_in_triangle(point, &triangle, &tri_heights)
                    {
                        raster[raster_idx] = Some(interpolated_height);
                    }
                }
            }
        }
    }

    raster
}

pub fn export_to_tiff(
    qm: &QuantizedMesh,
    filename: &Path,
    scale_shift: u16,
) -> Result<(), Box<dyn Error>> {
    let raster_size = get_raster_size(scale_shift);
    let raster = rasterize(qm, scale_shift);

    let mut tiff: TiffEncoder<File> = TiffEncoder::new(File::create(filename)?)?;
    let mut image = tiff
        .new_image::<Gray32Float>(raster_size.into(), raster_size.into())
        .unwrap();

    let geo_key_directory = vec![
    1, 1, 0, 2,          // Header: Version 1.1, 8 keys
    1024, 0, 1, 2,      // GTRasterTypeGeoKey 2 - geographic
    2048, 0, 1, 4979,   // GeodeticCRSGeoKey, value 4979 (WGS 84 3D)
    4096, 0, 1, 4979,   // VerticalGeoKey, value 4979 (WGS 84 3D)
    //0x1003, 0, 1, 9001,  // VerticalUnit, value 9001 (Meter)
   // 0x1005, 0, 1, 9102,  // HorizontalUnit, value 9102 (Degree)
   // 0x1004, 0, 1, 9102,  // GeographicUnit, value 9102 (Degree)
   // 0x2000, 0, 1, 0,     // Projection - set to 0 for geographic (or use appropriate EPSG)
    ];

    // Now define the bounding rectangle values
    let bounding_rectangle = &qm.bounding_rectangle; // Assuming this is your calculated bounding box
    let ll = bounding_rectangle.lower_left.to_degrees();
    let ur = bounding_rectangle.upper_right.to_degrees();

    let min_x = ll.lat();
    let min_y = ll.lon();
    let max_x = ur.lat();
    let max_y = ur.lon();
 
    // 0,0 is UL
    // crs axis order is lat,lon,alt
    // Define ModelTiepoint
    let tiepoints: Vec<f64> = vec![
        // 0.0,    // UL UV x
        // 0.0,    // UL UV y
        // 0.0,    // UL UV Z
        // *min_x, // UL world x
        // *max_y, // UL world y
        // 0.0,    // UL world z
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

    let image_data: Vec<f32> = raster
        .iter()
        .map(|&height| height.unwrap_or(f32::NAN)) // Replace None with NaN
        .collect();

    // Ensure we have the correct dimensions for the raster
    let expected_size = (raster_size as u32) * (raster_size as u32);
    if image_data.len() as u32 != expected_size {
        return Err(format!(
            "Raster size does not match expected dimensions: expected {}, found {}",
            expected_size,
            image_data.len()
        )
        .into());
    }
    //image.rows_per_strip(2).unwrap();

    // Write the image data as Gray32Float
    image.write_data(&image_data)?;
    Ok(())
}
