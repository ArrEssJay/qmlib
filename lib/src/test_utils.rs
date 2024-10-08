
// This module is only included during tests
#[cfg(test)]
pub mod test_data {
    use crate::geometry::{Ellipsoid, GeodeticPoint2};
    use crate::quantized_mesh_tile::{QuantizedMeshTile, TilingScheme, CRS};
    use crate::{geometry::CartesianPoint3, BoundingSphere, QuantizedMeshHeader};

    use crate::{EdgeIndices, QuantizedMesh, VertexData};

    // Test pattern generators
    // pyramid shape
    pub fn vd_test_pattern_pyramid() -> VertexData {
        let u = vec![0, 0, 32767, 32767, 16383];  // UV coordinates
        let v = vec![0, 32767, 32767, 0, 16383];
        let height = vec![0, 0, 0, 0, 32767];     // Heights

        let triangle_index = vec![
            [0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4],  // 4 triangles sharing the center vertex
        ];

        let vertex_count = u.len() as u32;
        let triangle_count = triangle_index.len() as u32;

        VertexData {
            vertex_count,
            u,
            v,
            height,
            triangle_count,
            triangle_index,
            edge_indices: EdgeIndices { west_vertex_count: 0, west_indices: vec![], south_vertex_count: 0, south_indices: vec![], east_vertex_count: 0, east_indices: vec![], north_vertex_count: 0, north_indices: vec![] },  // Edge indices can be added if needed
        }
    }

    // arrow pointing up/north
    pub fn vd_test_pattern_north() -> VertexData {
        let u = vec![0, 0, 16383, 32767, 32767];  // UV coordinates
        let v = vec![0, 32767, 32767, 32767, 0];
        let height = vec![0, 0, 32767, 0, 0];     // Heights

        let triangle_index = vec![
            [0, 1, 2], [2, 3, 4], [4, 0, 2],
        ];

        let vertex_count = u.len() as u32;
        let triangle_count = triangle_index.len() as u32;

        VertexData {
            vertex_count,
            u,
            v,
            height,
            triangle_count,
            triangle_index,
            edge_indices: EdgeIndices { west_vertex_count: 0, west_indices: vec![], south_vertex_count: 0, south_indices: vec![], east_vertex_count: 0, east_indices: vec![], north_vertex_count: 0, north_indices: vec![] },  // Edge indices can be added if needed
        }
    }

    fn vd_test_chessboard_mesh(grid_size: usize, high_value: u16) -> VertexData{
        let mut u = Vec::new();
        let mut v = Vec::new();
        let mut height = Vec::new();
        let mut triangle_index = Vec::new();
    
        // Step 1: Create grid points for u, v, and height
        let step = 32767 / grid_size as u16; // Normalizing UV coordinates to 32767
        for i in 0..=grid_size {
            for j in 0..=grid_size {
                u.push(i as u16 * step);
                v.push(j as u16 * step);
    
                // Alternate height for each square to create a checkerboard effect
                if (i + j) % 2 == 0 {
                    height.push(high_value); // High square
                } else {
                    height.push(0); // Low square
                }
            }
        }
    
        // Step 2: Create triangles for each square
        for i in 0..grid_size {
            for j in 0..grid_size {
                let top_left = i * (grid_size + 1) + j;
                let top_right = top_left + 1;
                let bottom_left = top_left + (grid_size + 1);
                let bottom_right = bottom_left + 1;
    
                // Create two triangles per square
                triangle_index.push([top_left, bottom_left, bottom_right]);  // First triangle
                triangle_index.push([top_left, bottom_right, top_right]);    // Second triangle
            }
        }
        let vertex_count = u.len() as u32;
        let triangle_count = triangle_index.len() as u32;

        VertexData {
            vertex_count,
            u,
            v,
            height,
            triangle_count,
            triangle_index,
            edge_indices: EdgeIndices { west_vertex_count: 0, west_indices: vec![], south_vertex_count: 0, south_indices: vec![], east_vertex_count: 0, east_indices: vec![], north_vertex_count: 0, north_indices: vec![] },  // Edge indices can be added if needed
        }

    }
    

    pub fn qm_header_test_1() -> QuantizedMeshHeader {
        QuantizedMeshHeader {
                center: CartesianPoint3::new(0.0,0.0,0.0),
                min_height: 0.0,
                max_height: 100.0,
                bounding_sphere: BoundingSphere { center: CartesianPoint3::new(0.0,0.0,0.0), radius: 50.0 },
                horizon_occlusion_point:CartesianPoint3::new(0.0,0.0,0.0),
            }
        }

    

    pub fn qm_test_pattern_north() -> QuantizedMesh {
        QuantizedMesh {
            header: qm_header_test_1(),
            vertex_data: vd_test_pattern_north(),
            extensions: vec![],
        }
    }

    pub fn qm_test_chessboard_mesh() -> QuantizedMesh {
        QuantizedMesh {
            header: qm_header_test_1(),
            vertex_data: vd_test_chessboard_mesh(32usize,16384u16),
            extensions: vec![],
        }
    }

   
    // override tile bounds crossing equator & prime meridian by 1 degree +/-
    pub fn qmt_test_north() ->QuantizedMeshTile {
        let mut qmt = QuantizedMeshTile::new(qm_test_pattern_north(), 15, 59489, 9692,Ellipsoid::wgs84(), TilingScheme::Tms, CRS::Epsg4326);
        qmt.bounding_rectangle.upper_right=GeodeticPoint2::from_degrees(1.0, 1.0);
        qmt.bounding_rectangle.lower_left=GeodeticPoint2::from_degrees(-1.0, -1.0);
        qmt
    } 

    pub fn qmt_test_chess() ->QuantizedMeshTile {
        let mut qmt = QuantizedMeshTile::new(qm_test_chessboard_mesh(), 15, 59489, 9692,Ellipsoid::wgs84(), TilingScheme::Tms, CRS::Epsg4326);
        qmt.bounding_rectangle.upper_right=GeodeticPoint2::from_degrees(1.0, 1.0);
        qmt.bounding_rectangle.lower_left=GeodeticPoint2::from_degrees(-1.0, -1.0);
        qmt
    } 
}
    

