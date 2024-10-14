use nalgebra::Point2;

use crate::geometry::Triangle;

pub enum InterpolationMethod {
    Edge,
    Barycentric,
}

pub fn interpolate_height(
    sample_point: &Point2<u16>,
    triangle: &Triangle<u16>,
    heights: &[f64; 3],
    method: &InterpolationMethod,
) -> Option<f32> {
    match method {
        InterpolationMethod::Edge => interpolate_height_edge(*sample_point, triangle, heights),
        InterpolationMethod::Barycentric => interpolate_height_barycentric(*sample_point, triangle, heights),
    }
}

pub fn interpolate_height_edge(
    sample_point: Point2<u16>,
    triangle: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    // Triangle vertices

    let mut v0x = triangle.vertices[0].x as i32;
    let mut v0y = triangle.vertices[0].y as i32;
    let v1x = triangle.vertices[1].x as i32;
    let v1y = triangle.vertices[1].y as i32;
    let mut v2x = triangle.vertices[2].x as i32;
    let mut v2y = triangle.vertices[2].y as i32;
    let px = sample_point.x as i32;
    let py = sample_point.y as i32;

     // Function to compute the edge function (signed area)
     let edge_function = |v0x: i32, v0y: i32, v1x: i32, v1y: i32, px: i32, py: i32| -> i32 {
        (px - v0x) * (v1y - v0y) - (py - v0y) * (v1x - v0x)
    };

    let is_ccw = |v0x: i32, v0y: i32, v1x: i32, v1y: i32, v2x: i32, v2y: i32| -> bool {
        // Using the determinant to check orientation
        ((v1x - v0x) * (v2y - v0y) - (v2x - v0x) * (v1y - v0y)) > 0
    };
    
    if !is_ccw(v0x, v0y, v1x, v1y, v2x, v2y) {
        let tx = v0x;
        let ty = v0y;
        v0x = v2x;
        v0y = v2y;
        v2x = tx;
        v2y = ty;
    }

    let w0 = edge_function(v0x, v0y, v1x, v1y, px, py);
    let w1 = edge_function(v1x, v1y, v2x, v2y, px, py);
    let w2 = edge_function(v2x, v2y, v0x, v0y, px, py);
    
    // Compute the full area of the triangle
    let area_full = edge_function(v0x, v0y, v1x, v1y, v2x, v2y);

    // If the area is zero, the triangle is degenerate
    if area_full == 0 {
        return None;
    }

    // Check if the point is inside the triangle (all weights must be non-negative)
    if w0 >= 0 && w1 >= 0 && w2 >= 0 {
        let numerator = w0  as f32 * heights[0] as f32 + w1 as f32 * heights[1] as f32 + w2 as f32 * heights[2] as f32;

        let interpolated_height = numerator / area_full as f32;

        Some(interpolated_height)
    } else {
        None // Point is outside the triangle
    }
}

pub fn interpolate_height_barycentric(
    point_usize: Point2<u16>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    // Cast triangle to f64 for bounds checking precision
    // Cast triangle vertices to f64 for precision
    let v0 = Point2::new(
        triangle_u16.vertices[0].x as f64,
        triangle_u16.vertices[0].y as f64,
    );
    let v1 = Point2::new(
        triangle_u16.vertices[1].x as f64,
        triangle_u16.vertices[1].y as f64,
    );
    let v2 = Point2::new(
        triangle_u16.vertices[2].x as f64,
        triangle_u16.vertices[2].y as f64,
    );

    let point: Point2<f64> = Point2::new(point_usize.x as f64, point_usize.y as f64);

    // Translate points so that v0 is the origin
    let p_prime = point - v0;
    let v1_prime = v1 - v0;
    let v2_prime = v2 - v0;

    // Create the 2x2 matrix A with v1_prime and v2_prime as columns
    let a = nalgebra::Matrix2::new(v1_prime.x, v2_prime.x, v1_prime.y, v2_prime.y);

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
            #[allow(clippy::cast_possible_truncation)]
            let interpolated_height =
                (lambda0 * heights[0] + lambda.x * heights[1] + lambda.y * heights[2]) as f32;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolate_height_edge_inside() {
        let triangle = Triangle {
            vertices: [Point2::new(0, 0), Point2::new(0, 10), Point2::new(10, 0)],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_edge(point, &triangle, &heights);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_edge_outside() {
        let triangle = Triangle {
            vertices: [Point2::new(10, 20), Point2::new(30, 40), Point2::new(50, 60)],
        };
        let heights = [100.0, 150.0, 200.0];
        let point = Point2::new(25, 30);
        let result = interpolate_height_edge(point, &triangle, &heights);
        assert_eq!(result, None);
    }

    
    #[test]
    fn test_interpolate_height_barycentric_inside() {
        let triangle = Triangle {
            vertices: [Point2::new(0, 0), Point2::new(10, 0), Point2::new(0, 10)],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_barycentric(point, &triangle, &heights);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_barycentric_outside() {
        let triangle = Triangle {
            vertices: [Point2::new(0, 0), Point2::new(10, 0), Point2::new(0, 10)],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_barycentric(point, &triangle, &heights);
        assert_eq!(result, None);
    }
}
