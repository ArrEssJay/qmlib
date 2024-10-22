
#[allow(unused_variables)]
#[allow(unused_attributes)]
#[allow(unused_imports)]
#[allow(dead_code)]
use std::cmp::{min, max};


use glam::{UVec2, UVec3};

pub struct ShaderUniforms {
    raster_dim_size: u32,
    height_min: f32,
    height_max: f32,
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
    for i in 0..vertices.len() {
        let v0 = vertices[indices[i][0] as usize];
        let v1 = vertices[indices[i][1] as usize];
        let v2 = vertices[indices[i][2] as usize];

        // Determine the bounding box of the triangle
        let min_x = min(min(v0.x, v1.x), v2.x);
        let min_y = min(min(v0.y, v1.y), v2.y);
        
        let max_x = max(max(v0.x, v1.x), v2.x);
        let max_y = max(max(v0.y, v1.y), v2.y);

        let raster_x_origin = min_x;

        // Invert y-axis by subtracting raster_point.y from raster_dim_size - 1
        let raster_y = uniforms.raster_dim_size - 1 - min_y;

        // Build 2D point as its needed for vector math in the triangle rasteriser
        let  mut raster_point: UVec2 = UVec2::new(raster_x_origin, raster_y);

        // NaN bit pattern
        //let f32_bits_nan = f32::to_bits(f32::NAN);

        for _ in 0..(max_y - min_y)
        // range is inclusive..exclusive
        {
            raster_point.x = raster_x_origin; // reset scanline x each row
            for _ in 0..(max_x - min_x) {
                // index in the flat raster
                let raster_idx = ((raster_y * uniforms.raster_dim_size) + raster_point.x) as usize;

                // Check if the raster cell is empty by reading the atomic value without locking
                //let cell = output[raster_idx as usize];

                if let Some(value) = shader_edge_interpolator(raster_point, [v0, v1, v2], uniforms) {
                    storage[raster_idx as usize] = value;
                }

                raster_point.x += 1;
            }
            raster_point.y -= 1;
        }
    }
    

}

// Not borrowing here as we are not modifying the input
pub fn shader_edge_interpolator(p: UVec2, v: [UVec3; 3], uniforms: &ShaderUniforms) -> Option<f32> {
    // Wrapping subtraction to avoid casting to signed integers
    // Taking x,y arguments as we have a mix of 2 & 3 dimensional vectors and this avoids casting
    let edge_function = |v0_x: u32, v0_y: u32, v1_x: u32, v1_y: u32, p_x: u32, p_y: u32| -> i32 {
        ((p_x.wrapping_sub(v0_x) as i32) * (v1_y.wrapping_sub(v0_y) as i32))
            - ((p_y.wrapping_sub(v0_y) as i32) * (v1_x.wrapping_sub(v0_x) as i32))
    };

    let is_ccw = |v0_x: u32, v0_y: u32, v1_x: u32, v1_y: u32, v2_x: u32, v2_y: u32| -> bool {
        // Using the determinant to check orientation
        ((v1_x.wrapping_sub(v0_x) as i32) * (v2_y.wrapping_sub(v0_y) as i32))
            - ((v2_x.wrapping_sub(v0_x) as i32) * (v1_y.wrapping_sub(v0_y) as i32))
            > 0
    };

    // check winding order is CCW and correct if necessary

    let (v0, v1, v2) = if is_ccw(
        v[0].x,
        v[0].y,
        v[1].x,
        v[1].y,
        v[2].x,
        v[2].y,
    ) {
        (v[0], v[1], v[2])
    } else {
        (v[2], v[1], v[0])
    };

    let w0 = edge_function(v0.x,v0.y, v1.x, v1.y, p.x,p.y);
    let w1 = edge_function(v1.x,v1.y,v2.x, v2.y, p.x,p.y);
    let w2 = edge_function(v2.x,v2.y,  v0.x,v0.y, p.x,p.y);

    // Compute the full area of the triangle
    let area_full = edge_function(v0.x,v0.y,v1.x,v1.y, v2.x, v2.y);

    // If the area is zero, the triangle is degenerate
    if area_full == 0 {
        return None;
    }

    // Check if the point is inside the triangle (all weights must be non-negative)
    // interpolate the z value
    if w0 >= 0 && w1 >= 0 && w2 >= 0 {
        let numerator = w0 as f32 * v0.z as f32 + w1 as f32 * v1.z as f32 + w2 as f32 * v2.z as f32;

        let interpolated_height = numerator / area_full as f32;

        // Normalize and map the height
        let normalized_height = interpolated_height / 32767.0;
        let mapped_height = uniforms.height_min + normalized_height * (uniforms.height_max - uniforms.height_min);

        Some(mapped_height)

    } else {
        None
    }
}
#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_shader_edge_interpolator_inside_triangle() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 100,
            height_min: 0.0,
            height_max: 1.0,
        };

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = shader_edge_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.5);
    }

    #[test]
    fn test_shader_edge_interpolator_outside_triangle() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 100,
            height_min: 0.0,
            height_max: 1.0,
        };

        let p = UVec2::new(3, 3);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = shader_edge_interpolator(p, v, &uniforms);
        assert!(result.is_none());
    }

    #[test]
    fn test_shader_edge_interpolator_on_edge() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 100,
            height_min: 0.0,
            height_max: 1.0,
        };

        let p = UVec2::new(1, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = shader_edge_interpolator(p, v, &uniforms);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.25);
    }

    #[test]
    fn test_shader_edge_interpolator_degenerate_triangle() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 100,
            height_min: 0.0,
            height_max: 1.0,
        };

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 0, 32767),
            UVec3::new(0, 0, 32767),
        ];

        let result = shader_edge_interpolator(p, v, &uniforms);
        assert!(result.is_none());
    }


    #[test]
    fn test_raster() {
        let uniforms = ShaderUniforms {
            raster_dim_size: 256,
            height_min: 0.0,
            height_max: 1.0,
        };

        let vertices = vec![
            UVec3::new(10, 10, 100),
            UVec3::new(20, 10, 200),
            UVec3::new(15, 20, 150),
        ];

        let indices = vec![
            [0, 1, 2],
        ];

        let mut storage = vec![0.0; (uniforms.raster_dim_size * uniforms.raster_dim_size) as usize];

        raster(&uniforms, &vertices, &indices, &mut storage);

        // Verify the output
        // This is a simple check to see if the storage buffer has been modified
        // You can add more specific checks based on the expected output
        assert!(storage.iter().any(|&value| value != 0.0));

        // Print the 2D array for visual inspection
        print_2d_array(&storage, uniforms.raster_dim_size as usize);
    }

    fn print_2d_array(storage: &[f32], raster_dim_size: usize) {
        for row in storage.chunks(raster_dim_size) {
            println!("{:?}", row);
        }
    }
}




