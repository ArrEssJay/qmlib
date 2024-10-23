use core::cmp::{max, min};

use spirv_std::glam::{UVec2, UVec3};

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct AABB {
    min: UVec2,
    max: UVec2,
}

#[derive(Clone, PartialEq)]
pub struct BVHNode {
    bounding_box: AABB,
    left: Option<Box<BVHNode>>,
    right: Option<Box<BVHNode>>,
    triangles: Vec<usize>,
}

#[derive(Clone, PartialEq)]
pub struct BVH {
    root: BVHNode,
}

#[derive(PartialEq)]
pub enum Axis {
    X,
    Y,
}

// Bounding Volume Hierachy (BVH) implementation
// Allows for quick intersection tests between a block of pixels and a list of triangles
impl BVH {
    pub fn build(vertices: &[UVec3], indices: &[[u32; 3]]) -> Self {
        let mut triangle_indices: Vec<usize> = (0..indices.len()).collect();
        let root = Self::build_recursive(vertices, indices, &mut triangle_indices);
        BVH { root }
    }

    fn build_recursive(
        vertices: &[UVec3],
        indices: &[[u32; 3]],
        triangle_indices: &mut [usize],
    ) -> BVHNode {
        if triangle_indices.is_empty() {
            return BVHNode {
                bounding_box: AABB {
                    min: UVec2::new(u32::MAX, u32::MAX),
                    max: UVec2::new(u32::MIN, u32::MIN),
                },
                left: None,
                right: None,
                triangles: vec![],
            };
        }

        if triangle_indices.len() == 1 {
            let triangle = &indices[triangle_indices[0]];
            return BVHNode {
                bounding_box: calculate_aabb(vertices, triangle),
                left: None,
                right: None,
                triangles: triangle_indices.to_vec(),
            };
        }

        let bounding_box = triangle_indices.iter().fold(
            AABB {
                min: UVec2::new(u32::MAX, u32::MAX),
                max: UVec2::new(u32::MIN, u32::MIN),
            },
            |mut aabb, &i| {
                let triangle = &indices[i];
                let tri_aabb = calculate_aabb(vertices, triangle);
                aabb.min = aabb.min.min(tri_aabb.min);
                aabb.max = aabb.max.max(tri_aabb.max);
                aabb
            },
        );

        let axis =
            if bounding_box.max.x - bounding_box.min.x > bounding_box.max.y - bounding_box.min.y {
                Axis::X
            } else {
                Axis::Y
            };

        triangle_indices.sort_by_key(|&i| {
            let triangle = &indices[i];
            let centroid = (vertices[triangle[0] as usize]
                + vertices[triangle[1] as usize]
                + vertices[triangle[2] as usize])
                / 3;
            if axis == Axis::X {
                centroid.x
            } else {
                centroid.y
            }
        });

        let mid = triangle_indices.len() / 2;
        let (left_indices, right_indices) = triangle_indices.split_at_mut(mid);

        BVHNode {
            bounding_box,
            left: Some(Box::new(Self::build_recursive(
                vertices,
                indices,
                left_indices,
            ))),
            right: Some(Box::new(Self::build_recursive(
                vertices,
                indices,
                right_indices,
            ))),
            triangles: vec![],
        }
    }
    pub fn find_intersecting_triangles(&self, x: u32, y: u32, size: u32) -> Vec<(AABB, usize)> {
        let block_aabb = AABB {
            min: UVec2::new(x, y),
            max: UVec2::new(x + size - 1, y + size - 1),
        };
        let mut result = vec![];
        self.traverse(&self.root, &block_aabb, &mut result);
        result
    }

    fn traverse(&self, node: &BVHNode, block_aabb: &AABB, result: &mut Vec<(AABB, usize)>) {
        if !self.aabb_intersects(&node.bounding_box, block_aabb) {
            return;
        }

        if node.triangles.is_empty() {
            if let Some(ref left) = node.left {
                self.traverse(left, block_aabb, result);
            }
            if let Some(ref right) = node.right {
                self.traverse(right, block_aabb, result);
            }
        } else {
            for &triangle_index in &node.triangles {
                result.push((node.bounding_box, triangle_index));
            }
        }
    }

    fn aabb_intersects(&self, a: &AABB, b: &AABB) -> bool {
        a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y
    }

    pub fn print_tree(&self) {
        self.print_node(&self.root, 0);
    }

    fn print_node(&self, node: &BVHNode, depth: usize) {
        let indent = "  ".repeat(depth);
        println!("{}Node at depth {}: {:?}", indent, depth, node.bounding_box);
        if let Some(ref left) = node.left {
            println!("{}  Left:", indent);
            self.print_node(left, depth + 1);
        }
        if let Some(ref right) = node.right {
            println!("{}  Right:", indent);
            self.print_node(right, depth + 1);
        }
        if !node.triangles.is_empty() {
            println!("{}  Triangles: {:?}", indent, node.triangles);
        }
    }
}

fn calculate_aabb(vertices: &[UVec3], indices: &[u32; 3]) -> AABB {
    let v0 = vertices[indices[0] as usize];
    let v1 = vertices[indices[1] as usize];
    let v2 = vertices[indices[2] as usize];
    AABB {
        min: UVec2::new(min(min(v0.x, v1.x), v2.x), min(min(v0.y, v1.y), v2.y)),
        max: UVec2::new(max(max(v0.x, v1.x), v2.x), max(max(v0.y, v1.y), v2.y)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_aabb_intersection() {
        let a = AABB {
            min: UVec2::new(0, 0),
            max: UVec2::new(5, 5),
        };
        let b = AABB {
            min: UVec2::new(3, 3),
            max: UVec2::new(8, 8),
        };
        let c = AABB {
            min: UVec2::new(6, 6),
            max: UVec2::new(10, 10),
        };

        let bvh = BVH {
            root: BVHNode {
                bounding_box: a,
                left: None,
                right: None,
                triangles: vec![],
            },
        };

        assert!(bvh.aabb_intersects(&a, &b));
        assert!(!bvh.aabb_intersects(&a, &c));
    }

    #[test]
    fn test_aabb_intersection_edge_touching() {
        let a = AABB {
            min: UVec2::new(0, 0),
            max: UVec2::new(5, 5),
        };
        let b = AABB {
            min: UVec2::new(5, 5),
            max: UVec2::new(10, 10),
        };

        let bvh = BVH {
            root: BVHNode {
                bounding_box: a,
                left: None,
                right: None,
                triangles: vec![],
            },
        };

        assert!(bvh.aabb_intersects(&a, &b));
    }
    #[test]
    fn test_aabb_intersection_no_overlap() {
        let a = AABB {
            min: UVec2::new(0, 0),
            max: UVec2::new(5, 5),
        };
        let b = AABB {
            min: UVec2::new(6, 6),
            max: UVec2::new(10, 10),
        };

        let bvh = BVH {
            root: BVHNode {
                bounding_box: a,
                left: None,
                right: None,
                triangles: vec![],
            },
        };

        assert!(!bvh.aabb_intersects(&a, &b));
    }

    #[test]
    fn test_aabb_intersection_contained() {
        let a = AABB {
            min: UVec2::new(0, 0),
            max: UVec2::new(10, 10),
        };
        let b = AABB {
            min: UVec2::new(3, 3),
            max: UVec2::new(7, 7),
        };

        let bvh = BVH {
            root: BVHNode {
                bounding_box: a,
                left: None,
                right: None,
                triangles: vec![],
            },
        };

        assert!(bvh.aabb_intersects(&a, &b));
    }

    #[test]
    fn test_aabb_intersection_identical() {
        let a = AABB {
            min: UVec2::new(0, 0),
            max: UVec2::new(5, 5),
        };

        let bvh = BVH {
            root: BVHNode {
                bounding_box: a,
                left: None,
                right: None,
                triangles: vec![],
            },
        };

        assert!(bvh.aabb_intersects(&a, &a));
    }

    #[test]
    fn test_empty_bvh() {
        let vertices = vec![];
        let indices = vec![];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(1, 1, 1);
        assert_eq!(intersecting_triangles.len(), 0);
    }

    #[test]
    fn test_build_recursive_single_triangle() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
        ];
        let indices = vec![[0, 1, 2]];
        let mut triangle_indices = vec![0];

        let node = BVH::build_recursive(&vertices, &indices, &mut triangle_indices);

        assert_eq!(node.triangles.len(), 1);
        assert_eq!(node.triangles[0], 0);
        assert!(node.left.is_none());
        assert!(node.right.is_none());
    }

    #[test]
    fn test_build_recursive_multiple_triangles() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
            UVec3::new(1, 1, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let mut triangle_indices = vec![0, 1];

        let node = BVH::build_recursive(&vertices, &indices, &mut triangle_indices);

        assert!(node.triangles.is_empty());
        assert!(node.left.is_some());
        assert!(node.right.is_some());
    }

    #[test]
    fn test_build_recursive_balanced_split() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 0),
            UVec3::new(0, 2, 0),
            UVec3::new(2, 2, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let mut triangle_indices = vec![0, 1];

        let node = BVH::build_recursive(&vertices, &indices, &mut triangle_indices);

        assert!(node.triangles.is_empty());
        assert!(node.left.is_some());
        assert!(node.right.is_some());
    }

    #[test]
    fn test_build_recursive_unbalanced_split() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(3, 0, 0),
            UVec3::new(0, 3, 0),
            UVec3::new(3, 3, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let mut triangle_indices = vec![0, 1];

        let node = BVH::build_recursive(&vertices, &indices, &mut triangle_indices);

        assert!(node.triangles.is_empty());
        assert!(node.left.is_some());
        assert!(node.right.is_some());
    }
    #[test]
    fn test_find_intersecting_triangles_single_point() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(4, 0, 0),
            UVec3::new(0, 4, 0),
            UVec3::new(4, 4, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(2, 2, 0);
        assert_eq!(intersecting_triangles.len(), 2);
    }
    #[test]
    fn test_find_intersecting_triangles() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 0),
            UVec3::new(0, 2, 0),
            UVec3::new(2, 2, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(1, 1, 1);
        assert_eq!(intersecting_triangles.len(), 2);
    }

    #[test]
    fn test_find_intersecting_triangles_no_intersection() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 0),
            UVec3::new(0, 2, 0),
            UVec3::new(2, 2, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(3, 3, 1);
        assert_eq!(intersecting_triangles.len(), 0);
    }

    #[test]
    fn test_find_intersecting_triangles_partial_intersection() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(2, 0, 0),
            UVec3::new(0, 2, 0),
            UVec3::new(2, 2, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(1, 1, 2);
        assert_eq!(intersecting_triangles.len(), 2);
    }

    #[test]
    fn test_find_intersecting_triangles_edge_touching() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(4, 0, 0),
            UVec3::new(0, 4, 0),
            UVec3::new(4, 4, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(2, 2, 1);
        assert_eq!(intersecting_triangles.len(), 2);
    }

    #[test]
    fn test_degenerate_triangles() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(0, 0, 0),
            UVec3::new(0, 0, 0),
        ];
        let indices = vec![[0, 1, 2]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(0, 0, 1);
        assert_eq!(intersecting_triangles.len(), 1);
    }

    #[test]
    fn test_bvh() {
        let mut vertices = vec![];
        let mut indices = vec![];
        for i in 1..500 {
            vertices.push(UVec3::new(i, i, 0));
            vertices.push(UVec3::new(i + 1, i, 0));
            vertices.push(UVec3::new(i, i + 1, 0));
            indices.push([(i - 1) * 3, (i - 1) * 3 + 1, (i - 1) * 3 + 2]);
        }

        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(10, 10, 1);
        assert_eq!(intersecting_triangles.len(), 2);
    }

    // Below are essentially the same tests as above, but with the additional check of the bounding box and
    // triangle parameters rather than simply their presence in the result
    #[test]
    fn test_find_intersecting_triangles_bbox() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
            UVec3::new(1, 1, 0),
            UVec3::new(2, 2, 0),
            UVec3::new(3, 2, 0),
            UVec3::new(2, 3, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3], [4, 5, 6]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(0, 0, 2);
        assert_eq!(intersecting_triangles.len(), 2);

        for (aabb, triangle_index) in intersecting_triangles {
            let triangle = &indices[triangle_index];
            let expected_aabb = calculate_aabb(&vertices, triangle);
            assert_eq!(aabb, expected_aabb);
        }
    }
    #[test]
    fn test_find_intersecting_triangles_partial_overlap_bbox() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
            UVec3::new(1, 1, 0),
            UVec3::new(2, 2, 0),
            UVec3::new(3, 2, 0),
            UVec3::new(2, 3, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3], [4, 5, 6]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(1, 1, 2);
        assert_eq!(intersecting_triangles.len(), 3);

        for (aabb, triangle_index) in intersecting_triangles {
            let triangle = &indices[triangle_index];
            let expected_aabb = calculate_aabb(&vertices, triangle);
            assert_eq!(aabb, expected_aabb);
        }
    }

    #[test]
    fn test_find_intersecting_triangles_no_overlap_bbox() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
            UVec3::new(1, 1, 0),
            UVec3::new(2, 2, 0),
            UVec3::new(3, 2, 0),
            UVec3::new(2, 3, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3], [4, 5, 6]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(10, 10, 2);
        assert_eq!(intersecting_triangles.len(), 0);
    }

    #[test]
    fn test_find_intersecting_triangles_edge_touching_bbox() {
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(1, 0, 0),
            UVec3::new(0, 1, 0),
            UVec3::new(1, 1, 0),
            UVec3::new(2, 2, 0),
            UVec3::new(3, 2, 0),
            UVec3::new(2, 3, 0),
        ];
        let indices = vec![[0, 1, 2], [1, 2, 3], [4, 5, 6]];
        let bvh = BVH::build(&vertices, &indices);

        let intersecting_triangles = bvh.find_intersecting_triangles(1, 1, 1);
        assert_eq!(intersecting_triangles.len(), 2);

        for (aabb, triangle_index) in intersecting_triangles {
            let triangle = &indices[triangle_index];
            let expected_aabb = calculate_aabb(&vertices, triangle);
            assert_eq!(aabb, expected_aabb);
        }
    }
}
