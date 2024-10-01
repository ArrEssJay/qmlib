use nalgebra::Point2;
use std::{error::Error, fs::File, path::Path};
use tiff::encoder::{colortype::Gray32Float, TiffEncoder};

use crate::{geometry::Triangle, QuantizedMesh, UV_MAX_U16};

pub fn get_raster_size(scale_shift: u16) -> u16 {
    (UV_MAX_U16 + 1) >> scale_shift
}

pub fn rasterize(qm: &QuantizedMesh, scale_shift: u16) -> Vec<Option<f32>> {
    let raster_size = get_raster_size(scale_shift);
    let mut raster: Vec<Option<f32>> = vec![None; raster_size as usize * raster_size as usize];
    let v = &qm.vertex_data;

    for triangle_index in 0..v.triangle_index.len() {
        let triangle_indices = &v.triangle_index[triangle_index];
        assert!(triangle_indices.iter().all(|&i| i < v.u.len()));

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

        let heights: [u16; 3] = [
            qm.vertex_data.height[triangle_indices[0] as usize],
            qm.vertex_data.height[triangle_indices[1] as usize],
            qm.vertex_data.height[triangle_indices[2] as usize],
        ];

        let triangle_bounds = triangle.bounding_rect();

        for px in (triangle_bounds.lower_left.x as u16).min(raster_size - 1)
            ..=(triangle_bounds.upper_right.x as u16).min(raster_size - 1)
        {
            for py in (triangle_bounds.lower_left.y as u16).min(raster_size - 1)
                ..=(triangle_bounds.upper_right.y as u16).min(raster_size - 1)
            {
                let point = Point2::new(px as f64, py as f64);

                if let Some(interpolated_height) =
                    barycentric_coords_and_interpolate_height(point, &triangle, &heights)
                {
                    let index = (py as usize * raster_size as usize) + px as usize;
                    raster[index] = Some(interpolated_height);
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

fn barycentric_coords_and_interpolate_height(
    point: Point2<f64>,
    t: &Triangle<f64>,
    heights: &[u16; 3],
) -> Option<f32> {
    let area = 0.5
        * ((t.vertices[1].x - t.vertices[0].x) * (t.vertices[2].y - t.vertices[0].y)
            - (t.vertices[2].x - t.vertices[0].x) * (t.vertices[1].y - t.vertices[0].y))
            .abs();

    let area_a = 0.5
        * ((t.vertices[1].x - point.x) * (t.vertices[2].y - point.y)
            - (t.vertices[2].x - point.x) * (t.vertices[1].y - point.y))
            .abs();

    let area_b = 0.5
        * ((t.vertices[2].x - point.x) * (t.vertices[0].y - point.y)
            - (t.vertices[0].x - point.x) * (t.vertices[2].y - point.y))
            .abs();

    let area_c = area - area_a - area_b;

    if area <= 0.0 {
        return None;
    }

    let a = area_a / area;
    let b = area_b / area;
    let c = area_c / area;

    if a >= 0.0 && b >= 0.0 && c >= 0.0 {
        let height = a as f32 * heights[0] as f32
            + b as f32 * heights[1] as f32
            + c as f32 * heights[2] as f32;
        Some(height)
    } else {
        None
    }
}
