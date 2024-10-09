use std::{fs, path::PathBuf};

use criterion::{criterion_group, criterion_main, Criterion};
use qmlib::{interpolator::{interpolate_height_barycentric, interpolate_height_lu, interpolate_height_parametric}, quantized_mesh_tile, tiff_writer::write_tiff};

fn benchmark(c: &mut Criterion) {
    let path: PathBuf = PathBuf::from(format!("{}/../test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));
    

    // Now a whole mesh

    let scale_shift: u16 = 4;
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();
    let output_path = path.with_extension("bench.tiff");

    c.bench_function("write_tiff_barycentric", |b| {
        b.iter(|| {
            write_tiff(&tile, &output_path, scale_shift, interpolate_height_barycentric).unwrap();
        })
    });

    c.bench_function("write_tiff_parametric", |b| {
        b.iter(|| {
            write_tiff(&tile, &output_path, scale_shift, interpolate_height_parametric).unwrap();
        })
    });

    c.bench_function("write_tiff_lu", |b| {
        b.iter(|| {
            write_tiff(&tile, &output_path, scale_shift, |point, triangle, heights| {
                interpolate_height_lu(point, triangle, heights, true) // Set bounds_check to true
            }).unwrap();
        })
    });

    // Clean up test output
    fs::remove_file(output_path).expect("Failed to clean up test output file");
}

criterion_group!(benches, benchmark);
criterion_main!(benches);