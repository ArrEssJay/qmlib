use std::path::PathBuf;

use criterion::{criterion_group, criterion_main, Criterion};
use qmlib::{interpolator::InterpolationMethod, quantized_mesh_tile, raster::rasterise};

fn benchmark(c: &mut Criterion) { 
    let path: PathBuf = PathBuf::from(format!("{}/../test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));
    

    // Now a whole mesh

    let raster_scale_factor: u16 = 4;
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();

    c.bench_function("rasterise_barycentric", |b| {
        b.iter(|| {
            rasterise(
                &tile,
                raster_scale_factor,
                &InterpolationMethod::Barycentric,
            )
        })
    });

    c.bench_function("write_tiff_edge", |b| {
        b.iter(|| {
            rasterise(
                &tile,
                raster_scale_factor,
                &InterpolationMethod::Edge,
            )
        })
    });

}

criterion_group!(benches, benchmark);
criterion_main!(benches);