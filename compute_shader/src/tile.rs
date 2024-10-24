use spirv_std::glam::{UVec2, UVec3, UVec4};

pub const GRID_CELL_SIZE: u32 = 8;

pub struct RasterParameters {
    pub raster_dim_size: u32,
    pub height_min: f32,
    pub height_max: f32,
}

#[allow(dead_code)] // Used in tests
impl RasterParameters {
    pub fn new(raster_dim_size: u32, height_min: f32, height_max: f32) -> Self {
        Self {
            raster_dim_size,
            height_min,
            height_max,
        }
    }
}

// Hack to make it easier to use UVec4 as an AABB
// We can't implement traits for types in other crates
// We can't use newtype pattern
// We can't even put it in a struct
// So we just define some functions that take a UVec4 as an AABB
// and return the values we need
// This is a bit ugly, but it works
pub type AABB = UVec4;
pub fn aabb_min_x(aabb: &AABB) -> u32 {
    aabb.x
}
pub fn aabb_min_y(aabb: &AABB) -> u32 {
    aabb.y
}
pub fn aabb_max_x(aabb: &AABB) -> u32 {
    aabb.z
}
pub fn aabb_max_y(aabb: &AABB) -> u32 {
    aabb.w
}
// Axis Aligned Bounding Box (AABB)
// #[derive(Clone, Copy, PartialEq)]
// pub struct AABB(UVec4);


// Just a convenience wrapper around UVec4
// impl AABB {
//     pub fn new(min_x: u32, min_y: u32, max_x: u32, max_y: u32) -> Self {
//         AABB(UVec4::new(min_x, min_y, max_x, max_y))
//     }

//     pub fn min_x(&self) -> u32 {
//         self.0.x
//     }

//     pub fn min_y(&self) -> u32 {
//         self.0.y
//     }

//     pub fn max_x(&self) -> u32 {
//         self.0.z
//     }

//     pub fn max_y(&self) -> u32 {
//         self.0.w
//     }

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

//... And a wrapper around UVec4::div()
// impl Div<u32> for AABB {
//     type Output = AABB;

//     fn div(self, rhs: u32) -> Self::Output {
//         AABB(self.0 / rhs)
//     }
// }

#[allow(dead_code)] // Used in tests
pub fn flat_cell_index(cell: UVec2, params: &RasterParameters) -> usize {
    (cell.y * (params.raster_dim_size / GRID_CELL_SIZE) + cell.x) as usize
}


// Only needed on host
// This doesn't really belong here but it avoids having to deal with
// importing from the shader crate
#[cfg(not(target_arch = "spirv"))]
pub fn assign_grid_cell_bounding_boxes(
    vertices: &[UVec3],
    indices: &[[u32; 3]], // vertex indices
) -> Vec<UVec4> {
 
    let mut bounding_boxes: Vec<UVec4> = Vec::new();

    for triangle_indices in indices.iter() {
        let aabb:UVec4 = calculate_triangle_aabb(vertices, triangle_indices);
        bounding_boxes.push(aabb);
    }

    bounding_boxes
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_new() {
        let aabb:AABB = AABB::new(1, 2, 3, 4);
        assert_eq!(aabb_min_x(&aabb), 1);
        assert_eq!(aabb_min_y(&aabb), 2);
        assert_eq!(aabb_max_x(&aabb), 3);
        assert_eq!(aabb_max_y(&aabb), 4);
    }

    #[test]
    fn test_aabb_div() {
        let aabb = AABB::new(8, 16, 24, 32);
        let result = aabb / 8;
        assert_eq!(aabb_min_x(&result), 1);
        assert_eq!(aabb_min_y(&result), 2);
        assert_eq!(aabb_max_x(&result), 3);
        assert_eq!(aabb_max_y(&result), 4);
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
        assert_eq!(aabb_min_x(&aabb), 0);
        assert_eq!(aabb_min_y(&aabb), 0);
        assert_eq!(aabb_max_x(&aabb), 10);
        assert_eq!(aabb_max_y(&aabb), 15);
    }

    #[test]
    fn test_assign_triangles_to_grid() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(10, 10, 0),
            UVec3::new(5, 15, 0),
        ];
        let indices= vec![[0, 1, 2]];

        let bounding_boxes = assign_grid_cell_bounding_boxes(&vertices, &indices);

        assert_eq!(bounding_boxes.len(), 1);
        let aabb = &bounding_boxes[0];
        assert_eq!(aabb_min_x(aabb), 0);
        assert_eq!(aabb_min_y(aabb), 0);
        assert_eq!(aabb_max_x(aabb), 10);
        assert_eq!(aabb_max_y(aabb), 15);
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
        assert_eq!(aabb_min_x(aabb1), 0);
        assert_eq!(aabb_min_y(aabb1), 0);
        assert_eq!(aabb_max_x(aabb1), 10);
        assert_eq!(aabb_max_y(aabb1), 15);

        let aabb2 = &bounding_boxes[1];
        assert_eq!(aabb_min_x(aabb2), 20);
        assert_eq!(aabb_min_y(aabb2), 20);
        assert_eq!(aabb_max_x(aabb2), 30);
        assert_eq!(aabb_max_y(aabb2), 30);
    }
    
}
