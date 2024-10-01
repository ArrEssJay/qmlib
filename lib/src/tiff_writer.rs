use nalgebra::{Matrix2, Point2};
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
                let raster_idx = (py as usize * raster_size as usize) + px as usize;
                if raster[raster_idx] == None {
                    if let Some(interpolated_height) =
                        barycentric_coords_and_interpolate_height(point, &triangle, &heights)
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

fn barycentric_coords_and_interpolate_height(
    point: Point2<f64>,
    triangle: &Triangle<f64>,
    heights: &[u16; 3], // Heights at the vertices
) -> Option<f32> {
    let v0 = triangle.vertices[0];
    let v1 = triangle.vertices[1];
    let v2 = triangle.vertices[2];

    // Translate points so that v0 is the origin
    let p_prime = point - v0;
    let v1_prime = v1 - v0;
    let v2_prime = v2 - v0;

    // Create the 2x2 matrix A with v1_prime and v2_prime as columns
    let a = Matrix2::new(v1_prime.x, v2_prime.x, v1_prime.y, v2_prime.y);

    // Check if the matrix is invertible (non-zero determinant)
    if let Some(a_inv) = a.try_inverse() {
        // Solve for lambda1 and lambda2
        let lambda = a_inv * p_prime;

        // Calculate lambda0
        let lambda0 = 1.0 - lambda.x - lambda.y;

        // Ensure that the barycentric coordinates are valid (inside the triangle)

        if lambda0 >= -f64::EPSILON && lambda.x >= -f64::EPSILON && lambda.y >= -f64::EPSILON {
            // Interpolate the height using the barycentric coordinates
            let interpolated_height = (lambda0 * heights[0] as f64
                + lambda.x * heights[1] as f64
                + lambda.y * heights[2] as f64) as f32;

            Some(interpolated_height)
        } else {
            // The point is outside the triangle
            None
        }
    } else {
        // The triangle is degenerate (area is zero), no solution
        None
    }
}
