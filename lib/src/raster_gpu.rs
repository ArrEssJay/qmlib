
use glam::UVec3;
use crate::quantized_mesh_tile;

use quantized_mesh_tile::QuantizedMeshTile;
use compute_shader_interface::{run_compute_shader, RasterParameters};



pub fn rasterise(
    qmt: &QuantizedMeshTile,
    scale_shift: u16,
    raster_dim_size: u32,

) {

    // Downsample X,Y vertices and zip into flattened vertex buffer
    let vertices: Vec<UVec3> = qmt.quantized_mesh.vertex_data.u.iter()
    .zip(qmt.quantized_mesh.vertex_data.v.iter())
    .zip(qmt.quantized_mesh.vertex_data.height.iter())
    .map(|((u, v), h)| UVec3::new((*u >> scale_shift).into(), (*v >> scale_shift).into(), (*h).into()))
    .collect();


    // Take every 3 sequential triangle indices and add to a buffer
    let indices: Vec<[u32; 3]> = qmt.quantized_mesh.vertex_data.triangle_index.iter()
    .map(|triangle| {
        let triangle: &[u32] = triangle.as_slice();
        [triangle[0], triangle[1], triangle[2]]
    })
    .collect();

    let params = RasterParameters {
        raster_dim_size: raster_dim_size,
        height_min: qmt.quantized_mesh.header.min_height,
        height_max:  qmt.quantized_mesh.header.max_height,
    };

    async_std::task::block_on(async {
        let result = run_compute_shader(
            &vertices,
            &indices,
            &params
        ).await;
        println!("Result length: {}", result.len());
    });


    
}
