use spirv_std::glam::{DVec2, DVec3, IVec2, UVec2, UVec3, Vec3Swizzles};

pub struct ShaderUniforms {
    raster_dim_size: u32,
    height_min: f64,
    height_max: f64,
}

#[allow(dead_code)]
// Used in tests
impl ShaderUniforms {
    fn new(raster_dim_size: u32, height_min: f64, height_max: f64) -> Self {
        Self {
            raster_dim_size,
            height_min,
            height_max,
        }
    }
}

// won't be run in this context
pub fn main() {
    
}

// compute shader logic, run on the CPU, for testing
pub fn raster(
     uniforms: &ShaderUniforms,
   vertices: &[UVec3],
     indices: &[[u32; 3]],
     storage: &mut [f32],
) {
    for i in 0..indices.len() {
        let v0 = vertices[indices[i][0] as usize];
        let v1 = vertices[indices[i][1] as usize];
        let v2 = vertices[indices[i][2] as usize];

        // Determine the bounding box of the triangle
        let ( min_x,  min_y,  max_x,  max_y): (u32,u32,u32,u32);
        #[cfg(not(target_arch = "spirv"))]
        {
            use std::cmp::{max, min};

            min_x = min(min(v0.x, v1.x), v2.x);
            min_y = min(min(v0.y, v1.y), v2.y);
            
            max_x = max(max(v0.x, v1.x), v2.x);
            max_y = max(max(v0.y, v1.y), v2.y);
        }
        #[cfg(target_arch = "spirv")]
        {
            use spirv_std::arch::{unsigned_max, unsigned_min};

            min_x = unsigned_min(unsigned_min(v0.x, v1.x), v2.x);
            min_y = unsigned_min(unsigned_min(v0.y, v1.y), v2.y);

            max_x = unsigned_max(unsigned_max(v0.x, v1.x), v2.x);
            max_y = unsigned_max(unsigned_max(v0.y, v1.y), v2.y);
        }

        // for dim_size=32768 max_x = 32767, max_y = 32767
        assert!(max_x > 0);
        assert!(max_y > 0);
        assert!(max_x < uniforms.raster_dim_size);
        assert!(max_y < uniforms.raster_dim_size);


        // Build 2D point with the location of the top left corner of the bounding box
        // Inverted y-axis. Note the decrementing y line counter
        let  mut raster_point: UVec2 = UVec2::new(min_x,  uniforms.raster_dim_size - 1 - min_y);

        // NaN bit pattern
        //let f32_bits_nan = f32::to_bits(f32::NAN);

        for _i in min_y..=max_y
        // range is inclusive..exclusive
        {
            raster_point.x = min_x; // reset scanline x each row
            for _j in min_x..=max_x {
                // index in the flat raster
                let raster_idx = ((raster_point.y * uniforms.raster_dim_size) + raster_point.x) as usize;

                // Check if the raster cell is empty by reading the atomic value without locking
                if let Some(value) = triangle_face_height_interpolator(raster_point, [v0, v1, v2], uniforms) {
                    
                    storage[raster_idx as usize] =value as f32;
                }
                raster_point.x += 1;
            }
            raster_point.y = raster_point.y.saturating_sub(1); //avoid underflow on min_y=0
        }
    }
}

// Is v2 inside the edge formed by v0 and v1
pub fn edge_function(v: [IVec2;3]) -> i32 {
    ((v[2].x  - v[0].x ) * (v[1].y  - v[0].y ) - (v[2].y  - v[0].y ) * (v[1].x  - v[0].x )) as i32
}

// Are the vertices in cw order
// Use the determinant to check orientation
pub fn is_cw(v: [IVec2; 3]) -> bool {
    (v[1].x -v[0].x) * (v[2].y -v[0].y) > (v[2].x - v[0].x) * (v[1].y - v[0].y)
}

// Calculate the edge weights for a point p
pub fn calculate_edge_weights(v: [IVec2; 3], p: IVec2) -> [i32;3] {
      [edge_function([v[0], v[1], p]),
     edge_function([v[1], v[2], p]),
    edge_function([v[2], v[0], p])]
   
}

// Is the point p inside the triangle formed by vertices v
// Uses orientation of 3 edges to determine if the point is inside the triangle
pub fn point_in_triangle(v: [UVec3; 3], p: UVec2) -> bool {
    // Get the xy components of the vertices
    let v_xy = v.map(|v| v.xy().as_ivec2());

    // Correct winding order to CCW if necessary
    let v_xy_ccw = if is_cw(v_xy) {
        [v_xy[2], v_xy[1], v_xy[0]]
    } else {
        v_xy
    };

    let w = calculate_edge_weights(v_xy_ccw, p.as_ivec2());

    // Compute the full area of the triangle
    let area_full = edge_function(v_xy_ccw); //4

    // If the area is zero, the triangle is degenerate
    if area_full == 0 {
        return false;
    }

    return w[0] >= 0 && w[1] >= 0 && w[2] >= 0;
}

// Calculate the barycentric weights for a point p
pub fn calculate_barycentric_weights(v: [DVec2; 3], p: DVec2) -> [f64;3] {

    let area_abc = (v[1] - v[0]).perp_dot(v[2] - v[0]).abs();
    let area_pbc = (v[1] - p).perp_dot(v[2] - p).abs();
    let area_pca = (v[2] - p).perp_dot(v[0] - p).abs();
    let area_pab = (v[0] - p).perp_dot(v[1] - p).abs();

    [ area_pbc / area_abc,
    area_pca / area_abc,
    area_pab / area_abc]

}

// Interpolate the height of a point p inside the triangle formed by vertices v
pub fn interpolate_barycentric(v: [DVec3; 3], p: DVec2,  uniforms: &ShaderUniforms) -> Option<f64> {
    
    let wb = calculate_barycentric_weights(v.map(|v| v.xy()), p);

    // Not checking weights as we have already decided that the point is inside the triangle
    let numerator = wb[0]* v[0].z + wb[1] * v[1].z + wb[2] * v[2].z; // 131068

    // Normalize and map the height. Unmapped range is 0-32767
    let normalized_height = numerator / 32767.;
    let mapped_height = uniforms.height_min + normalized_height * (uniforms.height_max - uniforms.height_min);
    Some(mapped_height)
}


pub fn triangle_face_height_interpolator(p: UVec2, v: [UVec3; 3], uniforms: &ShaderUniforms) -> Option<f64> {
   
    // Check if the point is inside the triangle (all weights must be non-negative)
    // using integer edge function weights. This avoids floating point precision issues.
    // Interpolate the z value using barycentric coordinates only if the point
    // is determined to be inside the triangle
    if point_in_triangle(v, p) {
        // Interpolate the z value using barycentric coordinates
        // work in double precision and reduce for output
        interpolate_barycentric(v.map(|v| v.as_dvec3()), p.as_dvec2(), uniforms)
    } else {
        None
    }
}

#[cfg(test)]

mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;


    #[test]
    fn test_is_cw() {
        // Clockwise points
        assert!(is_cw([IVec2::new(0, 0), IVec2::new(1, 0), IVec2::new(0, 1)]));
    
        // Counter-clockwise points
        assert!(!is_cw([IVec2::new(0, 0), IVec2::new(0, 1), IVec2::new(1, 0)]));
    
        // Collinear points (should not be considered CW)
        assert!(!is_cw([IVec2::new(0, 0), IVec2::new(1, 1), IVec2::new(2, 2)]));
    }
    
    #[test]
    fn test_edge_function_colinear() {
        // Point on the edge
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(4, 4);
        let v2 = IVec2::new(2, 2);
        let result = edge_function([v0, v1, v2]);
        assert_eq!(result, 0);
    }
    
    #[test]
    fn test_edge_function_left_of_edge() {
        // Point to the left of the edge
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(4, 4);
        let v2 = IVec2::new(1, 3);
        let result = edge_function([v0, v1, v2]);
        assert_eq!(result, -8);
    }
    
    #[test]
    fn test_edge_function_right_of_edge() {
        // Point to the right of the edge
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(4, 4);
        let v2 = IVec2::new(3, 1);
        let result = edge_function([v0, v1, v2]);
        assert_eq!(result, 8);
    }




    #[test]
    fn test_calculate_weights_inside_triangle() {
        // Note the CCW winding order
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(0, 3);
        let v2 = IVec2::new(3, 0);
        let p = IVec2::new(1, 1);

        let w = calculate_edge_weights([v0, v1, v2], p);
        assert_eq!(w[0], 3);
        assert_eq!(w[1], 3);
        assert_eq!(w[2], 3);
    }

    #[test]
    fn test_calculate_weights_outside_triangle() {
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(0, 2);
        let v2 = IVec2::new(2, 0);
        let p = IVec2::new(3, 3);

        let w = calculate_edge_weights([v0, v1, v2], p);
        assert_eq!(w[0], 6);
        assert_eq!(w[1], -8);
        assert_eq!(w[2], 6);
    }

    #[test]
    fn test_calculate_weights_on_edge() {
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(0, 2);
        let v2 = IVec2::new(2, 0);
        let p = IVec2::new(1, 1);

        let w = calculate_edge_weights([v0, v1, v2], p);
        assert_eq!(w[0], 2);
        assert_eq!(w[1], 0);
        assert_eq!(w[2], 2);
    }

    #[test]
    fn test_calculate_weights_at_vertex() {
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(0, 2);
        let v2 = IVec2::new(2, 0);
        let p = IVec2::new(0, 0);

        let w = calculate_edge_weights([v0, v1, v2], p);
        assert_eq!(w[0], 0);
        assert_eq!(w[1], 4);
        assert_eq!(w[2], 0);
    }

    #[test]
    fn test_calculate_weights_epsilon_outside() {
        //this point is outside the triangle but is determined to be inside due to floating point precision
        //edge method correctly determines it to be outside
        let v0 = IVec2::new(0, 0);
        let v1 = IVec2::new(0, 32767);
        let v2 = IVec2::new(1, 16384);
        let p = IVec2::new(1, 16383);

        let w = calculate_edge_weights([v0, v1, v2], p);
        assert_eq!(w[0], 32767);
        assert_eq!(w[1], 1);
        assert_eq!(w[2], -1);
    }

    #[test]
    fn test_barycentric_weights_inside_triangle() {
        let v0 = DVec2::new(0., 0.);
        let v1 = DVec2::new(0., 3.);
        let v2 = DVec2::new(3., 0.);
        let p = DVec2::new(1., 1.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0],1./3.);
        assert_eq!(w[1],1./3.);
        assert_eq!(w[2],1./3.);

    }

    #[test]
    fn test_barycentric_weights_on_edge() {
        let v0 = DVec2::new(0., 0.);
        let v1 = DVec2::new(0., 3.);
        let v2 = DVec2::new(3., 0.);
        let p = DVec2::new(1., 0.);
        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0],2./3.);
        assert_eq!(w[1],0.);
        assert_eq!(w[2],1./3.);
    }

    #[test]
    fn test_barycentric_weights_epsilon_outside() {
        //this point is outside the triangle but is determined to be inside due to floating point precision
        let v0 = DVec2::new(0., 0.);
        let v1 = DVec2::new(0., 32767.);
        let v2 = DVec2::new(1., 16384.);
        let p = DVec2::new(1., 16383.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_abs_diff_eq!(w[0],0., epsilon = 1e-4);
        assert_abs_diff_eq!(w[1],0., epsilon = 1e-4);
        assert_abs_diff_eq!(w[2],1., epsilon = 1e-5);
    }

    #[test]
    fn test_barycentric_weights_at_vertex() {
        let v0 = DVec2::new(0., 0.);
        let v1 = DVec2::new(2., 0.);
        let v2 = DVec2::new(0., 2.);
        let p = DVec2::new(0., 0.);

        let w = calculate_barycentric_weights([v0,v1,v2], p);
        assert_eq!(w[0],1.);
        assert_eq!(w[1],0.);
        assert_eq!(w[2],0.);
    }

    #[test]
    fn test_barycentric_weights_outside_triangle() {
        let v0 = DVec2::new(0., 0.);
        let v1 = DVec2::new(2., 0.);
        let v2 = DVec2::new(0., 2.);
        let p = DVec2::new(5., 5.);

        let w = calculate_barycentric_weights([v0,v1,v2], p);
        assert_eq!(w[0],4.);
        assert_eq!(w[1],2.5);
        assert_eq!(w[2],2.5);
    }

    #[test]
    fn test_shader_edge_interpolator_inside_triangle_zero_plane() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 3, 0),
            UVec3::new(3, 0, 0),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_inside_triangle_different_z() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 3, 32767),
            UVec3::new(3, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_abs_diff_eq!(result.unwrap(), 0.5, epsilon = 1e-5);
    }

    #[test]
    fn test_shader_edge_interpolator_on_edge_different_z() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 2, 32767),
            UVec3::new(2, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_abs_diff_eq!(result.unwrap(), 0.75, epsilon = 1e-5);
    }


    #[test]
    fn test_shader_edge_interpolator_outside_triangle() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(3, 3);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_none());
    }

    #[test]
    fn test_shader_edge_interpolator_at_vertex_different_z() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(0, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 2, 32767),
            UVec3::new(2, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_on_edge() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(1, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.5);
    }

    #[test]
    fn test_shader_edge_interpolator_at_vertex() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(0, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_degenerate_triangle() {
        let uniforms = ShaderUniforms::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 0, 32767),
            UVec3::new(0, 0, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &uniforms);
        assert!(result.is_none());
    }


    #[test]
    // The vertex positions are expected to be pre-scaled so as to be in the range of the raster size
    fn test_raster_plane() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 4,
            height_min: 0.,
            height_max: 100.,
        };

        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(0, 3, 0),
            UVec3::new(3, 3, 32767),
            UVec3::new(3, 0, 32767),
        ]; 

        let indices = vec![
            [0, 1, 2],
            [0, 2, 3],
        ];

        let mut storage = vec![-1.; (uniforms.raster_dim_size * uniforms.raster_dim_size) as usize];

        raster(&uniforms, &vertices, &indices, &mut storage);
       
        let expected_output = vec![
        0.0, 33.333336, 66.66667, 100.0,
        0.0, 33.333336, 66.66667, 100.0,
        0.0, 33.333336, 66.66667, 100.0,
        0.0, 33.333336, 66.66667, 100.0,
    ];
    assert_abs_diff_eq!(expected_output.as_slice(), storage.as_slice(), epsilon = 1e-5);

    }

}