use spirv_std::glam::{UVec2, UVec3, UVec4};

// This is the magic number that determines the size of the grid cells
// Its simply well suited to typical GPU thread group sizes
// and the size of the raster
// It should be a power of 2, and not greater than the raster size
// Or, the raster cannot be smaller than the grid cell size
pub const GRID_CELL_SIZE: u32 = 8;
// pub const NATIVE_RASTER_SIZE: u32 = 32768;
// pub const MAX_SCALE_SHIFT: u16 = 12; // 32768>>12 = 8 = 1 grid cell

pub struct RasterParameters {
    pub raster_dim_size: u32,
    pub height_min: f32,
    pub height_max: f32,
}

impl RasterParameters {
    pub fn new(raster_dim_size: u32, height_min: f32, height_max: f32) -> Self {
        Self {
            raster_dim_size,
            height_min,
            height_max,
        }
    }
}

// The spir-v compiler is very picky about how we extend/wrap
// types in external crates. traits on a type alias is the only way
// I've found to make this work
pub type AABB = UVec4;

pub trait AABBValues {
    fn min_x(&self) -> u32;
    fn min_y(&self) -> u32;
    fn max_x(&self) -> u32;
    fn max_y(&self) -> u32;
}

impl AABBValues for UVec4 {
    fn min_x(&self) -> u32 {
        self.x
    }

    fn min_y(&self) -> u32 {
        self.y
    }

    fn max_x(&self) -> u32 {
        self.z
    }

    fn max_y(&self) -> u32 {
        self.w
    }
}

// Only needed on host - requires std
#[cfg(not(target_arch = "spirv"))]

pub fn assign_grid_cell_bounding_boxes(
    vertices: &[UVec3],
    indices: &[[u32; 3]], // vertex indices
) -> Vec<UVec4> {
    let mut bounding_boxes: Vec<UVec4> = Vec::new();

    for triangle_indices in indices.iter() {
        let aabb: UVec4 = calculate_triangle_aabb(vertices, triangle_indices);
        bounding_boxes.push(aabb);
    }

    bounding_boxes
}

// Only needed on host - requires std
#[cfg(not(target_arch = "spirv"))]
pub fn calculate_triangle_aabb(vertices: &[UVec3], indices: &[u32; 3]) -> UVec4 {
    let v0 = vertices[indices[0] as usize];
    let v1 = vertices[indices[1] as usize];
    let v2 = vertices[indices[2] as usize];

    use core::cmp::{max, min};

    let min_x = min(min(v0.x, v1.x), v2.x);
    let min_y = min(min(v0.y, v1.y), v2.y);
    let max_x = max(max(v0.x, v1.x), v2.x);
    let max_y = max(max(v0.y, v1.y), v2.y);

    UVec4::new(min_x, min_y, max_x, max_y)
}

//  spir-v implementation of aabb - not needed at the moment
// #[cfg(target_arch = "spirv")]
// pub fn calculate_triangle_aabb(vertices: &UVec3, indices: &[u32; 3]) -> AABB {
//     let v0 = vertices[indices[0] as usize];
//     let v1 = vertices[indices[1] as usize];
//     let v2 = vertices[indices[2] as usize];

//     use spirv_std::arch::{unsigned_max, unsigned_min};

//     let min_x = unsigned_min(unsigned_min(v0.x, v1.x), v2.x);
//     let min_y = unsigned_min(unsigned_min(v0.y, v1.y), v2.y);
//     let max_x = unsigned_max(unsigned_max(v0.x, v1.x), v2.x);
//     let max_y = unsigned_max(unsigned_max(v0.y, v1.y), v2.y);

//     AABB(UVec4::new(min_x, min_y, max_x, max_y))
// }

//    }


#[allow(dead_code)] // Used in tests
pub fn flat_cell_index(cell: UVec2, params: &RasterParameters) -> usize {
    (cell.y * (params.raster_dim_size / GRID_CELL_SIZE) + cell.x) as usize
}



#[cfg(not(target_arch = "spirv"))]
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_new() {
        let aabb: AABB = AABB::new(1, 2, 3, 4);
        assert_eq!(aabb.min_x(), 1);
        assert_eq!(aabb.min_y(), 2);
        assert_eq!(aabb.max_x(), 3);
        assert_eq!(aabb.max_y(), 4);
    }

    #[test]
    fn test_aabb_div() {
        let aabb = AABB::new(8, 16, 24, 32);
        let result = aabb / 8;
        assert_eq!(result.min_x(), 1);
        assert_eq!(result.min_y(), 2);
        assert_eq!(result.max_x(), 3);
        assert_eq!(result.max_y(), 4);
    }

    #[test]
    fn test_calculate_triangle_aabb() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(10, 10, 0),
            UVec3::new(5, 15, 0),
        ];
        let indices = [0, 1, 2];
        let aabb = calculate_triangle_aabb(&vertices, &indices);
        assert_eq!(aabb.min_x(), 0);
        assert_eq!(aabb.min_y(), 0);
        assert_eq!(aabb.max_x(), 10);
        assert_eq!(aabb.max_y(), 15);
    }

    #[test]
    fn test_assign_triangles_to_grid() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(10, 10, 0),
            UVec3::new(5, 15, 0),
        ];
        let indices = vec![[0, 1, 2]];

        let bounding_boxes = assign_grid_cell_bounding_boxes(&vertices, &indices);

        assert_eq!(bounding_boxes.len(), 1);
        let aabb = &bounding_boxes[0];
        assert_eq!(aabb.min_x(), 0);
        assert_eq!(aabb.min_y(), 0);
        assert_eq!(aabb.max_x(), 10);
        assert_eq!(aabb.max_y(), 15);
    }

    #[test]
    fn test_flat_cell_index() {
        let params = RasterParameters::new(16, 0.0, 1.0);
        let cell = UVec2::new(1, 1);
        let index = flat_cell_index(cell, &params);
        assert_eq!(index, 3);
    }

    #[test]
    fn test_assign_grid_cell_bounding_boxes_multiple_triangles() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(10, 10, 0),
            UVec3::new(5, 15, 0),
            UVec3::new(20, 20, 0),
            UVec3::new(25, 25, 0),
            UVec3::new(30, 30, 0),
        ];
        let indices = vec![[0, 1, 2], [3, 4, 5]];

        let bounding_boxes = assign_grid_cell_bounding_boxes(&vertices, &indices);

        assert_eq!(bounding_boxes.len(), 2);

        let aabb1 = &bounding_boxes[0];
        assert_eq!(aabb1.min_x(), 0);
        assert_eq!(aabb1.min_y(), 0);
        assert_eq!(aabb1.max_x(), 10);
        assert_eq!(aabb1.max_y(), 15);

        let aabb2 = &bounding_boxes[1];
        assert_eq!(aabb2.min_x(), 20);
        assert_eq!(aabb2.min_y(), 20);
        assert_eq!(aabb2.max_x(), 30);
        assert_eq!(aabb2.max_y(), 30);
    }
}
