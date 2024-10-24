use std::fs::File;
use std::io::Read;

use spirv_std::glam::UVec3;
// use bytemuck::{cast_slice, Pod, Zeroable};

use wgpu::util::DeviceExt;
use wgpu::{Adapter, Features};

// #[repr(C)]
// #[derive(Copy, Clone, Pod,Zeroable)]

fn print_gpu_capabilities(adapter: &Adapter) {
    // Print adapter properties
    let adapter_info = adapter.get_info();
    println!("Adapter Info:");
    println!("  Name: {}", adapter_info.name);
    println!("  Vendor: {}", adapter_info.vendor);
    println!("  Device: {}", adapter_info.device);
    println!("  Type: {:?}", adapter_info.device_type);
    println!("  Backend: {:?}", adapter_info.backend);

    // Print supported features
    let features = adapter.features();
    println!("Supported Features:");
    for feature in wgpu::Features::all().iter() {
        if features.contains(feature) {
            println!("  {:?}", feature);
        }
    }

    // Print supported limits
    let limits = adapter.limits();
    println!("Supported Limits:");
    println!(
        "  Max Texture Dimension 1D: {}",
        limits.max_texture_dimension_1d
    );
    println!(
        "  Max Texture Dimension 2D: {}",
        limits.max_texture_dimension_2d
    );
    println!(
        "  Max Texture Dimension 3D: {}",
        limits.max_texture_dimension_3d
    );
    println!(
        "  Max Texture Array Layers: {}",
        limits.max_texture_array_layers
    );
    println!("  Max Bind Groups: {}", limits.max_bind_groups);
    println!(
        "  Max Dynamic Uniform Buffers Per Pipeline Layout: {}",
        limits.max_dynamic_uniform_buffers_per_pipeline_layout
    );
    println!(
        "  Max Dynamic Storage Buffers Per Pipeline Layout: {}",
        limits.max_dynamic_storage_buffers_per_pipeline_layout
    );
    println!(
        "  Max Sampled Textures Per Shader Stage: {}",
        limits.max_sampled_textures_per_shader_stage
    );
    println!(
        "  Max Samplers Per Shader Stage: {}",
        limits.max_samplers_per_shader_stage
    );
    println!(
        "  Max Storage Buffers Per Shader Stage: {}",
        limits.max_storage_buffers_per_shader_stage
    );
    println!(
        "  Max Storage Textures Per Shader Stage: {}",
        limits.max_storage_textures_per_shader_stage
    );
    println!(
        "  Max Uniform Buffers Per Shader Stage: {}",
        limits.max_uniform_buffers_per_shader_stage
    );
    println!(
        "  Max Uniform Buffer Binding Size: {}",
        limits.max_uniform_buffer_binding_size
    );
    println!(
        "  Max Storage Buffer Binding Size: {}",
        limits.max_storage_buffer_binding_size
    );
    println!("  Max Vertex Buffers: {}", limits.max_vertex_buffers);
    println!("  Max Vertex Attributes: {}", limits.max_vertex_attributes);
    println!(
        "  Max Vertex Buffer Array Stride: {}",
        limits.max_vertex_buffer_array_stride
    );
    println!(
        "  Max Inter-Stage Shader Components: {}",
        limits.max_inter_stage_shader_components
    );
    println!(
        "  Max Compute Workgroup Storage Size: {}",
        limits.max_compute_workgroup_storage_size
    );
    println!(
        "  Max Compute Invocations Per Workgroup: {}",
        limits.max_compute_invocations_per_workgroup
    );
    println!(
        "  Max Compute Workgroup Size X: {}",
        limits.max_compute_workgroup_size_x
    );
    println!(
        "  Max Compute Workgroup Size Y: {}",
        limits.max_compute_workgroup_size_y
    );
    println!(
        "  Max Compute Workgroup Size Z: {}",
        limits.max_compute_workgroup_size_z
    );
    println!(
        "  Max Compute Workgroups Per Dimension: {}",
        limits.max_compute_workgroups_per_dimension
    );
}

pub struct RasterParameters {
    raster_dim_size: u32,
    height_min: f32,
    height_max: f32,
}

pub async fn run_compute_shader(
    vertices: &[UVec3],
    indices: &[[u32; 3]],
    params: &RasterParameters,
) -> Vec<f32> {
    let output_raster_size_bytes = (params.raster_dim_size as u64 * params.raster_dim_size as u64)
        * std::mem::size_of::<f32>() as u64;
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

    print_gpu_capabilities(&adapter);

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
    let vertex_bytes: Vec<u8> = vertices
        .iter()
        .flat_map(|v| {
            v.to_array()
                .iter()
                .flat_map(|&x| x.to_ne_bytes())
                .collect::<Vec<_>>()
        })
        .collect();

    // Convert UVec3 to bytes
    // let index_bytes: Vec<u8> = cast_slice(&indices).to_vec();
    let index_bytes: Vec<u8> = indices
        .iter()
        .flat_map(|&triangle| {
            triangle
                .iter()
                .flat_map(|&x| x.to_ne_bytes())
                .collect::<Vec<_>>()
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

    let storage_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Storage Buffer"),
        size: output_raster_size_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    let readback_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Readback Buffer"),
        size: output_raster_size_bytes,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    // Convert each field to bytes
    let raster_dim_size_bytes = params.raster_dim_size.to_ne_bytes();
    let height_min_bytes = params.height_min.to_ne_bytes();
    let height_max_bytes = params.height_max.to_ne_bytes();

    // Concatenate the byte arrays
    let mut uniform_bytes = Vec::new();
    uniform_bytes.extend_from_slice(&raster_dim_size_bytes);
    uniform_bytes.extend_from_slice(&height_min_bytes);
    uniform_bytes.extend_from_slice(&height_max_bytes);

    let uniform_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Raster Parameters"),
        contents: &uniform_bytes,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
    });

    // Read the shader file at runtime

    let entry_point = "main_cs";
    //let compiled_shader_modules = maybe_watch();
    //let module = compiled_shader_modules.spv_module_for_entry_point(entry_point);

    // Read the SPIR-V file
    let path = "/Users/rowan/Projects-Code/qmlib/target/spirv-builder/spirv-unknown-vulkan1.1/release/deps/qmlib_compute_shader.spv";

    let mut file = File::open(path).expect("Failed to open SPIR-V file");
    let mut bytes = Vec::new();
    file.read_to_end(&mut bytes)
        .expect("Failed to read SPIR-V file");

    // Convert bytes to Vec<u32>
    let spirv: Vec<u32> = bytemuck::cast_slice(&bytes).to_vec();

    let label = "Shader Module";
    //let wgpu::ShaderModuleDescriptorSpirV { label, source } = module;
    let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some(label),
        source: wgpu::ShaderSource::SpirV(spirv.into()),
    });

    // bind group
    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: None,
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    has_dynamic_offset: false,
                    min_binding_size: None,
                    ty: wgpu::BufferBindingType::Uniform,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    has_dynamic_offset: false,
                    min_binding_size: None,
                    ty: wgpu::BufferBindingType::Storage { read_only: true }, // Set read_only to true
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 2,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    has_dynamic_offset: false,
                    min_binding_size: None,
                    ty: wgpu::BufferBindingType::Storage { read_only: true }, // Set read_only to true
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 3,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    has_dynamic_offset: false,
                    min_binding_size: None,
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                },
                count: None,
            },
        ],
    });

    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: uniform_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: vertex_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: index_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: storage_buffer.as_entire_binding(),
            },
        ],
        label: Some("Compute Bind Group"),
    });

    // pipeline
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("Compute Pipeline Layout"),
        bind_group_layouts: &[&bind_group_layout],
        push_constant_ranges: &[],
    });

    let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        compilation_options: Default::default(),
        cache: None,
        label: None,
        layout: Some(&pipeline_layout),
        module: &shader_module,
        entry_point: &entry_point,
    });

    let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
        label: Some("Command Encoder"),
    });

    // Calculate the number of workgroups needed
    let num_triangles = indices.len() as u32;
    let workgroup_size = 8 * 8; // 64 threads per workgroup
    let num_workgroups = (num_triangles + workgroup_size - 1) / workgroup_size;
    println!(
        "num_triangles: {}, workgroup_size: {}, num_workgroups:{}",
        num_triangles, workgroup_size, num_workgroups
    );

    // Set up the compute pass
    // Scope to ensure compute pass is dropped before the buffer is mapped
    {
        let mut compute_pass = encoder.begin_compute_pass(&Default::default());
        compute_pass.set_pipeline(&pipeline);
        compute_pass.set_bind_group(0, &bind_group, &[]);

        compute_pass.dispatch_workgroups(num_workgroups, 1, 1);
    }
    // Copy data from the storage buffer to the readback buffer
    encoder.copy_buffer_to_buffer(
        &storage_buffer,
        0,
        &readback_buffer,
        0,
        output_raster_size_bytes,
    );

    queue.submit(Some(encoder.finish()));

    // init data vec here
    let buffer_slice = readback_buffer.slice(..);

    buffer_slice.map_async(wgpu::MapMode::Read, |r| r.unwrap());

    //change this to add code to the
    //buffer_slice.map_async(wgpu::MapMode::Read, |r| r.unwrap());
    // NOTE(eddyb) `poll` should return only after the above callbacks fire
    // (see also https://github.com/gfx-rs/wgpu/pull/2698 for more details).
    device.poll(wgpu::Maintain::Wait);

    let data = buffer_slice.get_mapped_range();
    let result = data
        .chunks_exact(4)
        .map(|b| f32::from_ne_bytes(b.try_into().unwrap()))
        .collect::<Vec<_>>();
    drop(data);
    readback_buffer.unmap();

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[async_std::test]
    async fn test_run_compute_shader() {
        let params = RasterParameters {
            raster_dim_size: 64,
            height_min: 0.,
            height_max: 100.,
        };
        let vertices = vec![
            UVec3::new(0, 0, 0),
            UVec3::new(0, 63, 0),
            UVec3::new(63, 63, 32767),
            UVec3::new(63, 0, 32767),
        ];

        let indices = vec![[0, 1, 2], [0, 2, 3]];

        let result = run_compute_shader(&vertices, &indices, &params).await;

        // Add assertions to verify the result
        assert_eq!(
            result.len(),
            (params.raster_dim_size * params.raster_dim_size) as usize
        );
        let mut rows: Vec<Vec<f32>> = Vec::with_capacity(params.raster_dim_size as usize);

        for chunk in result.chunks(params.raster_dim_size as usize) {
            rows.push(chunk.to_vec());
        }

        for row in rows {
            println!("{:?}", row);
        }
    }
}
