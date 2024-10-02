use criterion::{criterion_group, criterion_main, Criterion};
use std::{fs, path::PathBuf};
use qmlib::{quantized_mesh_tile, tiff_writer};

fn benchmark(c: &mut Criterion) {
    let mut path: PathBuf = PathBuf::from(format!("{}/../test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));       
    let scale_shift: u16 = 4;
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();
    path.set_extension(".bench.tiff");

    c.bench_function("write_tiff", |b| b.iter(|| tiff_writer::write_tiff(&tile, &path, scale_shift).unwrap()));
    
    
    // Clean up test output
    fs::remove_file(path).expect("Failed to clean up test output file");
}

criterion_group!(benches, benchmark);
criterion_main!(benches);   