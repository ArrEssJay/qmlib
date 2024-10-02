use nalgebra::{Matrix2, Point2};

use crate::{geometry::Triangle, quantized_mesh_tile::QuantizedMeshTile, UV_MAX_U16};

pub fn raster_size_pixels(scale_shift: u16) -> u16 {
    (UV_MAX_U16 + 1) >> scale_shift
}

pub fn interpolate_triangle_face_height(
    point: Point2<f64>,
    triangle: &Triangle<f64>,
    heights: &[f64; 3],
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
        // Add a small epsilon threshold to account for floating-point precision errors.
        if lambda0 >= -f64::EPSILON && lambda.x >= -f64::EPSILON && lambda.y >= -f64::EPSILON {
            // Interpolate the height using the barycentric coordinates
            // Reduce precision to f32
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

pub fn rasterise(qmt: &QuantizedMeshTile, scale_shift: u16) -> Vec<Option<f32>> {
    let raster_size = raster_size_pixels(scale_shift);
    let mut raster: Vec<Option<f32>> = vec![None; raster_size as usize * raster_size as usize];
    let v = &qmt.quantized_mesh.vertex_data;

    let heights = &qmt.quantized_mesh.interpolated_height_vertices();
    for triangle_index in 0..v.triangle_index.len() {
        let triangle_indices = &v.triangle_index[triangle_index];

        let triangle = Triangle {
            vertices: [
                Point2::new(
                    (v.u[triangle_indices[0] as usize] >> scale_shift) as f64,
                    (v.v[triangle_indices[0] as usize] >> scale_shift) as f64,
                ),
                Point2::new(
                    (v.u[triangle_indices[1] as usize] >> scale_shift) as f64,
                    (v.v[triangle_indices[1] as usize] >> scale_shift) as f64,
                ),
                Point2::new(
                    (v.u[triangle_indices[2] as usize] >> scale_shift) as f64,
                    (v.v[triangle_indices[2] as usize] >> scale_shift) as f64,
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
                if raster[raster_idx] == None { // Only calculate once - triangles are non-overlapping 
                    if let Some(interpolated_height) =
                        interpolate_triangle_face_height(point, &triangle, &tri_heights)
                    {
                        raster[raster_idx] = Some(interpolated_height);
                    }
                }
            }
        }
    }

    raster
}