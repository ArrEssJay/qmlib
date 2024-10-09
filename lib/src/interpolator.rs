use nalgebra::{Matrix2, Matrix3, Point2, Vector2, Vector3};

use crate::geometry::Triangle;

pub type Interpolator = fn(Point2<usize>, &Triangle<u16>, &[f64; 3]) -> Option<f32>;

pub fn interpolate_height_lu(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    interpolate_height_lu_bounded(point_usize, triangle_u16, heights, true)
}

pub fn interpolate_height_lu_bounded(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
    bounds_check: bool,
) -> Option<f32> {
    // Cast triangle vertices to f64 for precision
    let v0 = Point2::new(triangle_u16.vertices[0].x as f64, triangle_u16.vertices[0].y as f64);
    let v1 = Point2::new(triangle_u16.vertices[1].x as f64, triangle_u16.vertices[1].y as f64);
    let v2 = Point2::new(triangle_u16.vertices[2].x as f64, triangle_u16.vertices[2].y as f64);

    // Create the matrix and vector for the system of equations
    let a = Matrix3::new(
        v0.x, v0.y, 1.0,
        v1.x, v1.y, 1.0,
        v2.x, v2.y, 1.0,
    );
    let b = Vector3::new(heights[0], heights[1], heights[2]);

    // Solve for the plane coefficients using QR decomposition
    if let Some(coefficients) = a.lu().solve(&b) {
        let a = coefficients[0];
        let b = coefficients[1];
        let c = coefficients[2];

        // Cast point to f64
        let point = Point2::new(point_usize.x as f64, point_usize.y as f64);

        // Calculate the height at the given point
        let z = a * point.x + b * point.y + c;

        if bounds_check {
            // Solve for s and t using the inverse of the matrix formed by v1 - v0 and v2 - v0
            let v0v1 = v1 - v0;
            let v0v2 = v2 - v0;
            let v0p = point - v0;

            let a = Matrix2::new(v0v1.x, v0v2.x, v0v1.y, v0v2.y);
            if let Some(inv_a) = a.try_inverse() {
                let st = inv_a * Vector2::new(v0p.x, v0p.y);
                let s = st.x;
                let t = st.y;

                // Check if the point is inside the triangle
                if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
                     Some(z as f32)
                } else {
                    // The point is outside the triangle
                    None
                }
            } else {
                // The matrix is singular, no solution
                None
            }
        } else {
            Some(z as f32)
        }
    }
    else {
        None
    }
}

pub fn interpolate_height_qr(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    interpolate_height_qr_bounded(point_usize, triangle_u16, heights, true)
}

pub fn interpolate_height_qr_bounded(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
    bounds_check: bool,
) -> Option<f32> {
    // Cast triangle vertices to f64 for precision
    let v0 = Point2::new(triangle_u16.vertices[0].x as f64, triangle_u16.vertices[0].y as f64);
    let v1 = Point2::new(triangle_u16.vertices[1].x as f64, triangle_u16.vertices[1].y as f64);
    let v2 = Point2::new(triangle_u16.vertices[2].x as f64, triangle_u16.vertices[2].y as f64);

    // Create the matrix and vector for the system of equations
    let a = Matrix3::new(
        v0.x, v0.y, 1.0,
        v1.x, v1.y, 1.0,
        v2.x, v2.y, 1.0,
    );
    let b = Vector3::new(heights[0], heights[1], heights[2]);

    // Solve for the plane coefficients using QR decomposition
    if let Some(coefficients) = a.qr().solve(&b) {
        let a = coefficients[0];
        let b = coefficients[1];
        let c = coefficients[2];

        // Cast point to f64
        let point = Point2::new(point_usize.x as f64, point_usize.y as f64);

        // Calculate the height at the given point
        let z = a * point.x + b * point.y + c;

        if bounds_check {
            // Solve for s and t using the inverse of the matrix formed by v1 - v0 and v2 - v0
            let v0v1 = v1 - v0;
            let v0v2 = v2 - v0;
            let v0p = point - v0;

            let a = Matrix2::new(v0v1.x, v0v2.x, v0v1.y, v0v2.y);
            if let Some(inv_a) = a.try_inverse() {
                let st = inv_a * Vector2::new(v0p.x, v0p.y);
                let s = st.x;
                let t = st.y;

                // Check if the point is inside the triangle
                if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
                     Some(z as f32)
                } else {
                    // The point is outside the triangle
                    None
                }
            } else {
                // The matrix is singular, no solution
                None
            }
        } else {
            Some(z as f32)
        }
    }
    else {
        None
    }
}

pub fn interpolate_height_parametric(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    // Cast triangle vertices to f64 for precision
    let v0 = Point2::new(triangle_u16.vertices[0].x as f64, triangle_u16.vertices[0].y as f64);
    let v1 = Point2::new(triangle_u16.vertices[1].x as f64, triangle_u16.vertices[1].y as f64);
    let v2 = Point2::new(triangle_u16.vertices[2].x as f64, triangle_u16.vertices[2].y as f64);

    // Cast point to f64
    let point = Point2::new(point_usize.x as f64, point_usize.y as f64);

    // Define vectors v and w
    let v = Vector2::new(v1.x - v0.x, v1.y - v0.y);
    let w = Vector2::new(v2.x - v0.x, v2.y - v0.y);

    // Define the point relative to v0
    let p = Vector2::new(point.x - v0.x, point.y - v0.y);

    // Solve for s and t using the inverse of the matrix formed by v and w
    let matrix = Matrix2::new(v.x, w.x, v.y, w.y);
    if let Some(params) = matrix.try_inverse().map(|inv| inv * p) {
        let s = params.x;
        let t = params.y;

        // Check if the point is inside the triangle
        if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
            // Interpolate the height
            let z = (1.0 - s - t) * heights[0] + s * heights[1] + t * heights[2];
            return Some(z as f32);
        }
    }
    None
}



pub fn interpolate_height_barycentric(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {
    // Cast triangle to f64 for bounds checking precision
    // Cast triangle vertices to f64 for precision
    let v0 = Point2::new(triangle_u16.vertices[0].x as f64, triangle_u16.vertices[0].y as f64);
    let v1 = Point2::new(triangle_u16.vertices[1].x as f64, triangle_u16.vertices[1].y as f64);
    let v2 = Point2::new(triangle_u16.vertices[2].x as f64, triangle_u16.vertices[2].y as f64);

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
            let interpolated_height = (lambda0 * heights[0]
                + lambda.x * heights[1] 
                + lambda.y * heights[2] ) as f32;

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
    fn test_interpolate_height_qr_inside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_qr_bounded(point, &triangle, &heights, true);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_qr_outside_bounded() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_qr_bounded(point, &triangle, &heights, true);
        assert_eq!(result, None);
    }

    #[test]
    fn test_interpolate_height_qr_outside_unbounded() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_qr_bounded(point, &triangle, &heights, false);
        assert_eq!(result, Some(45.0)); // note its outside but still calculcated
    }

    #[test]
    fn test_interpolate_height_parametric_inside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_parametric(point, &triangle, &heights);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_parametric_outside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_parametric(point, &triangle, &heights);
        assert_eq!(result, None);
    }

    #[test]
    fn test_interpolate_height_lu_inside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_lu_bounded(point, &triangle, &heights, true);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_lu_outside_bounded() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_lu_bounded(point, &triangle, &heights, true);
        assert_eq!(result, None);
    }

    #[test]
    fn test_interpolate_height_lu_outside_unbounded() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_lu_bounded(point, &triangle, &heights, false);
        assert_eq!(result, Some(45.0)); // note its outside but still calculcated
    }

    #[test]
    fn test_interpolate_height_barycentric_inside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(5, 5);
        let result = interpolate_height_barycentric(point, &triangle, &heights);
        assert_eq!(result, Some(15.0));
    }

    #[test]
    fn test_interpolate_height_barycentric_outside() {
        let triangle = Triangle {
            vertices: [
                Point2::new(0, 0),
                Point2::new(10, 0),
                Point2::new(0, 10),
            ],
        };
        let heights = [0.0, 10.0, 20.0];
        let point = Point2::new(15, 15);
        let result = interpolate_height_barycentric(point, &triangle, &heights);
        assert_eq!(result, None);
    }
}
