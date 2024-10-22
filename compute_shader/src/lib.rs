#![cfg_attr(target_arch = "spirv", no_std)]
// HACK(eddyb) can't easily see warnings otherwise from `spirv-builder` builds.
#![deny(warnings)]


use glam::{UVec2, UVec3};
use spirv_std::{
    arch::{unsigned_max, unsigned_min},
    glam, spirv,
};


pub struct ShaderUniforms {
    raster_dim_size: u32,
    height_min: f32,
    height_max: f32,
}

// Workgroup size is 8x8x1 (x,y,z)
#[spirv(compute(threads(256)))]
pub fn main_cs(
    #[spirv(uniform, descriptor_set = 0, binding = 0)] uniforms: &ShaderUniforms,
    #[spirv(storage_buffer, descriptor_set = 0, binding = 1)] vertices: &[UVec3],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 2)] indices: &[[u32; 3]],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 3)] output: &mut [f32],
) {
    for i in 0..vertices.len() {
        let v0 = vertices[indices[i][0] as usize];
        let v1 = vertices[indices[i][1] as usize];
        let v2 = vertices[indices[i][2] as usize];

        // Determine the bounding box of the triangle
        let min_x = unsigned_min(unsigned_min(v0.x, v1.x), v2.x);
        let min_y = unsigned_min(unsigned_min(v0.y, v1.y), v2.y);

        let max_x = unsigned_max(unsigned_max(v0.x, v1.x), v2.x);
        let max_y = unsigned_max(unsigned_max(v0.y, v1.y), v2.y);

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
                    output[raster_idx as usize] = value;
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
