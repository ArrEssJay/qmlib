use std::fs::File;
use std::io::Read;

use spirv_std::glam::UVec3;
// use bytemuck::{cast_slice, Pod, Zeroable};

use wgpu::util::DeviceExt;
use wgpu::Features;

// #[repr(C)]
// #[derive(Copy, Clone, Pod,Zeroable)]

pub struct ShaderPushConstants {
    pub raster_dim_size: u32,
    pub height_min: f32,
    pub height_max: f32,
}



pub async fn run_compute_shader(
    raster_dim_size: u32,
    vertices: &[UVec3],
    indices: &[[u32;3]],
    height_range: &[f32;2],
) -> Vec<u32> {
    let output_raster_size_bytes =
        (raster_dim_size as u64 * raster_dim_size as u64) * std::mem::size_of::<f32>() as u64;
    // device
    let backends = wgpu::util::backend_bits_from_env().unwrap_or(wgpu::Backends::PRIMARY);
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends,
        dx12_shader_compiler: wgpu::util::dx12_shader_compiler_from_env().unwrap_or_default(),
        ..Default::default()
    });
    let adapter = wgpu::util::initialize_adapter_from_env_or_default(&instance, None)
        .await
        .expect("Failed to find an appropriate adapter");

    let (device, queue) = adapter
        .request_device(
            &wgpu::DeviceDescriptor {
                label: None,
                required_features: Features::empty(),
                required_limits: wgpu::Limits::default(),
                memory_hints: wgpu::MemoryHints::Performance,
            },
            None,
        )
        .await
        .expect("Failed to create device");
    drop(instance);
    drop(adapter);

    // Convert UVec3 to bytes
    //let vertex_bytes: Vec<u8> = cast_slice(&vertices).to_vec();
    let vertex_bytes: Vec<u8> = vertices.iter()
    .flat_map(|v| {
        v.to_array().iter().flat_map(|&x| x.to_ne_bytes()).collect::<Vec<_>>()
    })
    .collect();

    // Convert UVec3 to bytes
   // let index_bytes: Vec<u8> = cast_slice(&indices).to_vec();
   let index_bytes: Vec<u8> = indices.iter()
   .flat_map(|&triangle| {
       triangle.iter().flat_map(|&x| x.to_ne_bytes()).collect::<Vec<_>>()
   })
   .collect();

    //buffers - map to byte arrays
    let vertex_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Vertex Buffer"),
        contents: vertex_bytes.as_slice(),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let index_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Index Buffer"),
        contents: index_bytes.as_slice(),
        usage: wgpu::BufferUsages::STORAGE,
    });

    let output_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Output Buffer"),
        size: output_raster_size_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    

    // Read the shader file at runtime

    let entry_point = "main_cs";
    //let compiled_shader_modules = maybe_watch();
    //let module = compiled_shader_modules.spv_module_for_entry_point(entry_point);
    
    // Read the SPIR-V file
    let path = "/Users/rowan/Projects-Code/qmlib/target/spirv-builder/spirv-unknown-vulkan1.1/release/deps/qmlib_compute_shader.spv";



     let mut file = File::open(path).expect("Failed to open SPIR-V file");
     let mut bytes = Vec::new();
     file.read_to_end(&mut bytes).expect("Failed to read SPIR-V file");
 
     // Convert bytes to Vec<u32>
     let spirv: Vec<u32> = bytemuck::cast_slice(&bytes).to_vec();
 
    let label ="Shader Module";
    //let wgpu::ShaderModuleDescriptorSpirV { label, source } = module;
    let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::SpirV(spirv.into()),
    });


    // Define the push constant range
    let push_constant_range = wgpu::PushConstantRange {
        stages: wgpu::ShaderStages::COMPUTE,
        range: 0..12, // 2 * f32 (4 bytes each) + 1 * u32 (4 bytes) = 12 bytes
    };

      // bind group
      let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: None,
        entries: &[wgpu::BindGroupLayoutEntry {
            binding: 0,
            count: None,
            visibility: wgpu::ShaderStages::COMPUTE,
            ty: wgpu::BindingType::Buffer {
                has_dynamic_offset: false,
                min_binding_size: None,
                ty: wgpu::BufferBindingType::Storage { read_only: false },
            },
        }],
    });

    // pipeline
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("Compute Pipeline Layout"),
        bind_group_layouts: &[&bind_group_layout],
        push_constant_ranges: &[push_constant_range],
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        compilation_options: Default::default(),
        cache: None,
        label: None,
        layout: Some(&pipeline_layout),
        module: &shader_module,
        entry_point: &entry_point,
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: vertex_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: index_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: output_buffer.as_entire_binding(),
            },
        ],
        label: Some("Compute Bind Group"),
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Compute Encoder"),
    });

    let mut compute_pass = encoder.begin_compute_pass(&Default::default());
    compute_pass.set_pipeline(&pipeline);
    compute_pass.set_bind_group(0, &bind_group, &[]);

    
    let push_constants = ShaderPushConstants {
        raster_dim_size: raster_dim_size,
        height_min: height_range[0],
        height_max: height_range[1],
    };

    // Convert each field to bytes
    let raster_dim_size_bytes = push_constants.raster_dim_size.to_ne_bytes();
    let height_min_bytes = push_constants.height_min.to_ne_bytes();
    let height_max_bytes = push_constants.height_max.to_ne_bytes();

    // Concatenate the byte arrays
    let mut push_constant_bytes = Vec::new();
    push_constant_bytes.extend_from_slice(&raster_dim_size_bytes);
    push_constant_bytes.extend_from_slice(&height_min_bytes);
    push_constant_bytes.extend_from_slice(&height_max_bytes);

   // compute_pass.set_push_constants(0, bytemuck::bytes_of(&push_constants));
    compute_pass.set_push_constants(0, &push_constant_bytes);


    compute_pass.dispatch_workgroups((raster_dim_size + 7) / 8, (raster_dim_size + 7) / 8, 1);

    queue.submit(Some(encoder.finish()));

    let buffer_slice = output_buffer.slice(..);
    buffer_slice.map_async(wgpu::MapMode::Read, |r| r.unwrap());
    // NOTE(eddyb) `poll` should return only after the above callbacks fire
    // (see also https://github.com/gfx-rs/wgpu/pull/2698 for more details).
    device.poll(wgpu::Maintain::Wait);

    let data = buffer_slice.get_mapped_range();
    //let result: Vec<u32> = cast_slice(&data).to_vec();

    let result: Vec<u32> = data.chunks_exact(4)
    .map(|chunk| u32::from_ne_bytes(chunk.try_into().expect("slice with incorrect length")))
    .collect();

    output_buffer.unmap();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[async_std::test]
    async fn test_run_compute_shader() {
        let raster_dim_size = 256;
        let height_range =[100.0f32,200.0f32];
        let vertices = vec![UVec3::new(0, 0,100), UVec3::new(0, 255,200), UVec3::new(255, 0,0)];
        let indices = [[0, 1, 2]];

        let result = run_compute_shader(raster_dim_size, &vertices, &indices, &height_range).await;

        // Add assertions to verify the result
        assert_eq!(result.len(), (raster_dim_size * raster_dim_size) as usize);
    }
}