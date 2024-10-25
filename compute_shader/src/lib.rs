#![cfg_attr(target_arch = "spirv", no_std)]
mod tile;

use spirv_std::{
    glam::{IVec2, UVec2, UVec3, UVec4, Vec2, Vec3, Vec3Swizzles},
    spirv,
};
use tile::{AABBValues, RasterParameters, AABB, GRID_CELL_SIZE};

#[allow(unused_imports)] // Spir-v compiler will complain if we don't
use spirv_std::num_traits::Float;

// If testing on the host we need these to dispatch the cells to be rasterised
#[cfg(not(target_arch = "spirv"))]
use tile::assign_grid_cell_bounding_boxes;

// Workgroup size is 8x8x1 (x,y,z)
// For now we are computing the bvh on the CPU and passing it to the shader
// Ideally we want to compute the bvh on the GPU as well, but that's for another day
// Workgroup size is 8x8x1 (x,y,z)
#[spirv(compute(threads(8, 8, 1)))]
pub fn main_cs(
    #[spirv(global_invocation_id)] global_id: UVec3,
    #[spirv(uniform, descriptor_set = 0, binding = 0)] params: &RasterParameters,
    #[spirv(storage_buffer, descriptor_set = 0, binding = 1)] vertices: &[UVec3],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 2)] indices: &[[u32; 3]],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 3)] bounding_boxes: &[UVec4],
    #[spirv(storage_buffer, descriptor_set = 0, binding = 4)] storage: &mut [f32],
) {
    let cell_index = global_id.xy();
    rasterise_cell(
        params,
        vertices,
        indices,
        bounding_boxes,
        storage,
        cell_index,
    );
}

// Build a raster using a scanline approach over each triangle's bounding box
// Only needed on host
#[cfg(not(target_arch = "spirv"))]
pub fn build_raster_triangle_scanline(
    params: &RasterParameters,
    vertices: &[UVec3],
    indices: &[[u32; 3]],
    storage: &mut [f32],
) {
    for index in 0..indices.len() {
        rasterise_triangle(params, vertices, indices, storage, index)
    }
}

// Build a raster by partitioning the raster into blocks, and rasterising the portion
// of the triangles contained within each block
// This first pre-computes triangle bounding boxes for the grid cell to use for
// intersection tests
// Only needed on host
#[cfg(not(target_arch = "spirv"))]
pub fn build_raster_block_scanline(
    params: &RasterParameters,
    vertices: &[UVec3],
    indices: &[[u32; 3]],
    storage: &mut [f32],
) {
    let bounding_boxes = assign_grid_cell_bounding_boxes(vertices, indices);
    let grid_size = params.raster_dim_size / GRID_CELL_SIZE;
    for y in 0..grid_size {
        for x in 0..grid_size {
            rasterise_cell(
                params,
                vertices,
                indices,
                &bounding_boxes,
                storage,
                UVec2::new(x, y),
            );
        }
    }
}

fn rasterise_cell(
    params: &RasterParameters,
    vertices: &[UVec3],
    indices: &[[u32; 3]],
    bounding_boxes: &[AABB],
    storage: &mut [f32],
    cell_index: UVec2,
) {
    // Calculate the base raster index for the cell
    let raster_index = cell_index * GRID_CELL_SIZE;

    for y in (raster_index.y..raster_index.y + GRID_CELL_SIZE).rev() {
        for x in raster_index.x..raster_index.x + GRID_CELL_SIZE {
            let pixel = UVec2::new(x, y);

            // Iterate over the bounding boxes and check for intersection
            for i in 0..bounding_boxes.len() {
                let aabb: AABB = bounding_boxes[i];
                let grid_aabb: AABB = aabb / GRID_CELL_SIZE;

                // Check if the bounding box intersects the current cell
                if grid_aabb.min_x() <= cell_index.x
                    && grid_aabb.max_x() >= cell_index.x
                    && grid_aabb.min_y() <= cell_index.y
                    && grid_aabb.max_y() >= cell_index.y
                {
                    let triangle_indices = indices[i];
                    let triangle = [
                        vertices[triangle_indices[0] as usize],
                        vertices[triangle_indices[1] as usize],
                        vertices[triangle_indices[2] as usize],
                    ];

                    // Edge function to determine if the point is inside the triangle
                    if point_in_triangle(triangle, pixel) {
                        if let Some(value) =
                            triangle_face_height_interpolator(pixel, triangle, params)
                        {
                            let raster_idx = (y * params.raster_dim_size + x) as usize;
                            // Only this block will write to these elements of the storage buffer so no need
                            // to care about atomicity
                            storage[raster_idx] = value;
                            break; // It's a keeper
                        }
                    }
                }
            }
        }
    }
}

// Only needed on host
#[cfg(not(target_arch = "spirv"))]
pub fn rasterise_triangle(
    params: &RasterParameters,
    vertices: &[UVec3],
    indices: &[[u32; 3]],
    storage: &mut [f32],
    index: usize,
) {
    use tile::calculate_triangle_aabb;

    let v0 = vertices[indices[index][0] as usize];
    let v1 = vertices[indices[index][1] as usize];
    let v2 = vertices[indices[index][2] as usize];

    let aabb: AABB = calculate_triangle_aabb(vertices, &indices[index]);

    // invert the y axis line order
    for y in (aabb.min_y()..=aabb.max_y()).rev() {
        for x in aabb.min_x()..=aabb.max_x() {
            // index in the flat raster
            let raster_idx = ((y * params.raster_dim_size) + x) as usize;

            // Check if the raster cell is empty by reading the atomic value without locking
            if let Some(value) =
                triangle_face_height_interpolator(UVec2::new(x, y), [v0, v1, v2], params)
            {
                // this invites a race condition as multiple threads can write to the same cell though in
                // theory they should be writing the same value
                storage[raster_idx] = value;
            }
        }
    }
}

// Is v2 inside the edge formed by v0 and v1
pub fn edge_function(v: [IVec2; 3]) -> i32 {
    ((v[2].x - v[0].x) * (v[1].y - v[0].y) - (v[2].y - v[0].y) * (v[1].x - v[0].x)) as i32
}

// Are the vertices in cw order
// Use the determinant to check orientation
pub fn is_cw(v: [IVec2; 3]) -> bool {
    (v[1].x - v[0].x) * (v[2].y - v[0].y) > (v[2].x - v[0].x) * (v[1].y - v[0].y)
}

// Calculate the edge weights for a point p
pub fn calculate_edge_weights(v: [IVec2; 3], p: IVec2) -> [i32; 3] {
    [
        edge_function([v[0], v[1], p]),
        edge_function([v[1], v[2], p]),
        edge_function([v[2], v[0], p]),
    ]
}

// Is the point p inside the triangle formed by vertices v
// Uses orientation of 3 edges to determine if the point is inside the triangle
pub fn point_in_triangle(v: [UVec3; 3], p: UVec2) -> bool {
    // Get the xy components of the vertices
    //let v_xy = v.map(|v| v.xy().as_ivec2());
    //RJ - cannot cast between pointer types -- spirv
    let v_xy: [IVec2; 3] = [
        v[0].xy().as_ivec2(),
        v[1].xy().as_ivec2(),
        v[2].xy().as_ivec2(),
    ];

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
pub fn calculate_barycentric_weights(v: [Vec2; 3], p: Vec2) -> [f32; 3] {
    //area of the sub-triangles formed by the vertices and point p
    let area_abc = f32::abs((v[1] - v[0]).perp_dot(v[2] - v[0]));
    let area_pbc = f32::abs((v[1] - p).perp_dot(v[2] - p));
    let area_pca = f32::abs((v[2] - p).perp_dot(v[0] - p));
    let area_pab = f32::abs((v[0] - p).perp_dot(v[1] - p));
    [
        area_pbc / area_abc,
        area_pca / area_abc,
        area_pab / area_abc,
    ]
}

// Interpolate the height of a point p inside the triangle formed by vertices v
pub fn interpolate_barycentric(v: [Vec3; 3], p: Vec2, params: &RasterParameters) -> Option<f32> {
    // RJ - error casting pointers spirv
    //let wb = calculate_barycentric_weights(v.map(|v| v.xy()), p);
    let v_xy: [Vec2; 3] = [v[0].xy(), v[1].xy(), v[2].xy()];
    let wb = calculate_barycentric_weights(v_xy, p);

    // Not checking weights as we have already decided that the point is inside the triangle
    let numerator = wb[0] * v[0].z + wb[1] * v[1].z + wb[2] * v[2].z; // 131068

    // Normalize and map the height. Unmapped range is 0-32767
    let normalized_height = numerator / 32767.;
    let mapped_height =
        params.height_min + normalized_height * (params.height_max - params.height_min);
    Some(mapped_height)
}

// Check if the point is inside the triangle (all weights must be non-negative)
// using integer edge function weights. This avoids floating point precision issues.
// Interpolate the z value using barycentric coordinates only if the point
// is determined to be inside the triangle
pub fn triangle_face_height_interpolator(
    p: UVec2,
    v: [UVec3; 3],
    params: &RasterParameters,
) -> Option<f32> {
    if point_in_triangle(v, p) {
        // Interpolate the z value using barycentric coordinates
        // Ideally work in double precision and reduce for output
        // As of now, f64 seems to be broken in rust-gpu

        // -- Error casting pointers in spirv compiler
        // interpolate_barycentric(v.map(|v| v.as_vec3()), p.as_vec2(), params)
        let v_vec3: [Vec3; 3] = [v[0].as_vec3(), v[1].as_vec3(), v[2].as_vec3()];
        interpolate_barycentric(v_vec3, p.as_vec2(), params)
    } else {
        None
    }
}

#[cfg(test)]

mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use paste::paste;

    // For the plane tests below this is simply a linear gradient from 0 to 100
    fn generate_expected_output_row(dim: usize) -> Vec<f32> {
        (0..dim)
            .map(|x| 100.0 * (x as f32 / (dim - 1) as f32))
            .collect()
    }

    macro_rules! generate_plane_raster_test {
        ($raster_fn:ident, $fn_prefix:ident, $dim:expr) => {
            paste! {
                #[test]
                fn [<$fn_prefix _ $dim>]() {
                    let params = RasterParameters {
                        raster_dim_size: $dim,
                        height_min: 0.,
                        height_max: 100.,
                    };

                    let vertices = vec![
                        UVec3::new(0, 0, 0),
                        UVec3::new(0, $dim as u32 - 1, 0),
                        UVec3::new($dim as u32 - 1, $dim as u32 - 1, 32767),
                        UVec3::new($dim as u32 - 1, 0, 32767),
                    ];

                    let indices = vec![[0, 2, 1], [0, 3, 2]];

                    let mut storage = vec![-1.; (params.raster_dim_size * params.raster_dim_size) as usize];

                    $raster_fn(&params, &vertices, &indices, &mut storage);

                    let expected_output_row = generate_expected_output_row($dim);

                    for row in 0..params.raster_dim_size {
                        let start = (row * params.raster_dim_size) as usize;
                        let end = start + params.raster_dim_size as usize;
                        let actual_row = &storage[start..end];
                        assert_abs_diff_eq!(
                            expected_output_row.as_slice(),
                            actual_row,
                            epsilon = 1e-4
                        );
                    }
                }
            }
        };
    }

    // Generate the test functions for given dimension and each rasteriser
    macro_rules! generate_plane_raster_block_scanline_tests {
        ($($dim:expr),*) => {
            $(
                generate_plane_raster_test!(build_raster_block_scanline, test_raster_block_scanline, $dim);
            )*
        };
    }

    macro_rules! generate_plane_triangle_rasteriser_tests {
        ($($dim:expr),*) => {
            $(
                generate_plane_raster_test!(build_raster_triangle_scanline, test_raster_triangle_scanline, $dim);
            )*
        };
    }



    generate_plane_triangle_rasteriser_tests!(8, 16, 32, 64, 128, 256, 512, 1024);
    //generate_plane_triangle_rasteriser_tests!(2048, 4096, 8192, 16384, 32768); // large
    generate_plane_raster_block_scanline_tests!(8, 16, 32, 64, 128, 256, 512, 1024);
    //generate_plane_raster_block_scanline_tests!(2048, 4096, 8192, 16384, 32768); //large


    #[test]
    fn test_is_cw() {
        // Clockwise winding
        assert!(is_cw([
            IVec2::new(0, 0),
            IVec2::new(1, 0),
            IVec2::new(0, 1)
        ]));

        // Counter-clockwise winding
        assert!(!is_cw([
            IVec2::new(0, 0),
            IVec2::new(0, 1),
            IVec2::new(1, 0)
        ]));

        // Collinear (should not be considered CW)
        assert!(!is_cw([
            IVec2::new(0, 0),
            IVec2::new(1, 1),
            IVec2::new(2, 2)
        ]));
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
        let v0 = Vec2::new(0., 0.);
        let v1 = Vec2::new(0., 3.);
        let v2 = Vec2::new(3., 0.);
        let p = Vec2::new(1., 1.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0], 1. / 3.);
        assert_eq!(w[1], 1. / 3.);
        assert_eq!(w[2], 1. / 3.);
    }

    #[test]
    fn test_barycentric_weights_on_edge() {
        let v0 = Vec2::new(0., 0.);
        let v1 = Vec2::new(0., 3.);
        let v2 = Vec2::new(3., 0.);
        let p = Vec2::new(1., 0.);
        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0], 2. / 3.);
        assert_eq!(w[1], 0.);
        assert_eq!(w[2], 1. / 3.);
    }

    #[test]
    fn test_barycentric_weights_epsilon_outside() {
        //this point is outside the triangle but is determined to be inside due to floating point precision
        let v0 = Vec2::new(0., 0.);
        let v1 = Vec2::new(0., 32767.);
        let v2 = Vec2::new(1., 16384.);
        let p = Vec2::new(1., 16383.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_abs_diff_eq!(w[0], 0., epsilon = 1e-4);
        assert_abs_diff_eq!(w[1], 0., epsilon = 1e-4);
        assert_abs_diff_eq!(w[2], 1., epsilon = 1e-5);
    }

    #[test]
    fn test_barycentric_weights_at_vertex() {
        let v0 = Vec2::new(0., 0.);
        let v1 = Vec2::new(2., 0.);
        let v2 = Vec2::new(0., 2.);
        let p = Vec2::new(0., 0.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0], 1.);
        assert_eq!(w[1], 0.);
        assert_eq!(w[2], 0.);
    }

    #[test]
    fn test_barycentric_weights_outside_triangle() {
        let v0 = Vec2::new(0., 0.);
        let v1 = Vec2::new(2., 0.);
        let v2 = Vec2::new(0., 2.);
        let p = Vec2::new(5., 5.);

        let w = calculate_barycentric_weights([v0, v1, v2], p);
        assert_eq!(w[0], 4.);
        assert_eq!(w[1], 2.5);
        assert_eq!(w[2], 2.5);
    }

    #[test]
    fn test_shader_edge_interpolator_inside_triangle_zero_plane() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 3, 0),
            UVec3::new(3, 0, 0),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_inside_triangle_different_z() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 3, 32767),
            UVec3::new(3, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_abs_diff_eq!(result.unwrap(), 0.5, epsilon = 1e-5);
    }

    #[test]
    fn test_shader_edge_interpolator_on_edge_different_z() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 2, 32767),
            UVec3::new(2, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_abs_diff_eq!(result.unwrap(), 0.75, epsilon = 1e-5);
    }

    #[test]
    fn test_shader_edge_interpolator_outside_triangle() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(3, 3);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_none());
    }

    #[test]
    fn test_shader_edge_interpolator_at_vertex_different_z() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(0, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 2, 32767),
            UVec3::new(2, 0, 16383),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_on_edge() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(1, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.5);
    }

    #[test]
    fn test_shader_edge_interpolator_at_vertex() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(0, 0);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 32767),
            UVec3::new(0, 2, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_some());
        assert_eq!(result.unwrap(), 0.);
    }

    #[test]
    fn test_shader_edge_interpolator_degenerate_triangle() {
        let params = RasterParameters::new(100, 0., 1.);

        let p = UVec2::new(1, 1);
        let v = [
            UVec3::new(0, 0, 0),
            UVec3::new(0, 0, 32767),
            UVec3::new(0, 0, 32767),
        ];

        let result = triangle_face_height_interpolator(p, v, &params);
        assert!(result.is_none());
    }

    // Both methods should produce an -identical- raster
    #[test]
    fn test_raster_methods_consistency() {
        let params = RasterParameters {
            raster_dim_size: 1024,
            height_min: 0.,
            height_max: 100.,
        };
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(0, 1023, 0),
            UVec3::new(1023, 1023, 32767),
            UVec3::new(1023, 0, 32767),
        ];

        let indices = vec![[0, 1, 2], [0, 2, 3]];

        let mut storage_triangle =
            vec![-1.; (params.raster_dim_size * params.raster_dim_size) as usize];
        let mut storage_bvh = vec![-1.; (params.raster_dim_size * params.raster_dim_size) as usize];

        // Build raster using triangle scanning method
        build_raster_triangle_scanline(&params, &vertices, &indices, &mut storage_triangle);

        // Build raster using block scanning method
        build_raster_block_scanline(&params, &vertices, &indices, &mut storage_bvh);

        // Compare the outputs
        assert_abs_diff_eq!(
            storage_triangle.as_slice(),
            storage_bvh.as_slice(),
            epsilon = 1e-4
        );
    }
}
