#![cfg_attr(target_arch = "spirv", no_std)]

use spirv_std::spirv;
use glam::{Vec2, Vec3, Vec4};

#[spirv(compute(threads(8, 8, 1)))]
pub fn main_cs(
    #[spirv(global_invocation_id)] id: Vec3,
    #[spirv(storage_buffer, descriptor_set = 0, binding = 0)] vertices: &[Vec4],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 1)] indices: &[u32],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 2)] output: &mut [u32],
) {
    let width = 800; // Assuming a fixed width for simplicity
    let height = 600; // Assuming a fixed height for simplicity

    let x = id.x as u32;
    let y = id.y as u32;

    if x >= width || y >= height {
        return;
    }

    let pixel_index = y * width + x;

    // Iterate over triangles
    for i in (0..indices.len()).step_by(3) {
        let v0 = vertices[indices[i] as usize];
        let v1 = vertices[indices[i + 1] as usize];
        let v2 = vertices[indices[i + 2] as usize];

        // Barycentric coordinates
        let p = Vec2::new(x as f32 + 0.5, y as f32 + 0.5);
        let v0p = Vec2::new(v0.x, v0.y) - p;
        let v1p = Vec2::new(v1.x, v1.y) - p;
        let v2p = Vec2::new(v2.x, v2.y) - p;

        let area = (v1.x - v0.x) * (v2.y - v0.y) - (v2.x - v0.x) * (v1.y - v0.y);
        let w0 = (v1p.x * v2p.y - v1p.y * v2p.x) / area;
        let w1 = (v2p.x * v0p.y - v2p.y * v0p.x) / area;
        let w2 = 1.0 - w0 - w1;

        if w0 >= 0.0 && w1 >= 0.0 && w2 >= 0.0 {
            // Interpolate z-value
            let z = w0 * v0.z + w1 * v1.z + w2 * v2.z;
            output[pixel_index as usize] = z.to_bits();
        }
    }
}