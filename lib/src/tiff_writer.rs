use nalgebra::Point2;
use std::{error::Error, fs::File, path::Path};
use tiff::encoder::{colortype::Gray32Float, TiffEncoder};

use crate::{geometry::Triangle, raster::interpolate_height_in_triangle, QuantizedMesh, UV_MAX_U16};

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

    let mut encoder: TiffEncoder<File> = TiffEncoder::new(File::create(filename)?)?;

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

    // Write the image data as Gray32Float
    encoder.write_image::<Gray32Float>(raster_size.into(), raster_size.into(), &image_data)?;

    Ok(())
}