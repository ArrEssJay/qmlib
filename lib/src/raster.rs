use nalgebra::{Matrix2, Point2};

use crate::geometry::Triangle;

pub fn interpolate_height_in_triangle(
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