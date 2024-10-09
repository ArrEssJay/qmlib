use std::{fs, path::PathBuf};

use criterion::{criterion_group, criterion_main, Criterion};
use qmlib::{interpolator::{ interpolate_height_barycentric, interpolate_height_edge, Interpolator}, quantized_mesh_tile, tiff_writer::write_tiff};
use qmlib::interpolator::InterpolationStrategy;

fn benchmark(c: &mut Criterion) { 
    let path: PathBuf = PathBuf::from(format!("{}/../test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));
    

    // Now a whole mesh

    let scale_shift: u16 = 4;
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();
    let output_path = path.with_extension("bench.tiff");

     // Define the plane-based interpolator and solver
    //  let plane_interpolator: PlaneInterpolator = |point, triangle, _heights, plane| {
    //     interpolate_height_plane(point, triangle, plane)
    // };
    // let lu_plane_solver: PlaneSolver = solve_plane_coefficients_lu;
    // let qr_plane_solver: PlaneSolver = solve_plane_coefficients_qr;

    let barycentric_interpolator: Interpolator = interpolate_height_barycentric;
    let edge_interpolator: Interpolator = interpolate_height_edge;

    c.bench_function("write_tiff_barycentric", |b| {
        b.iter(|| {
            write_tiff(
                &tile,
                &output_path,
                scale_shift,
                InterpolationStrategy::Simple(barycentric_interpolator),
            ).unwrap();
        })
    });

    c.bench_function("write_tiff_edge", |b| {
        b.iter(|| {
            write_tiff(
                &tile,
                &output_path,
                scale_shift,
                InterpolationStrategy::Simple(edge_interpolator),
            ).unwrap();
        })
    });

    // c.bench_function("write_tiff_qr_bounded", |b| {
    //     b.iter(|| {
    //         write_tiff(
    //             &tile,
    //             &output_path,
    //             scale_shift,
    //             InterpolationStrategy::PlaneBased {
    //                 plane_interpolator,
    //                 plane_solver:qr_plane_solver,
    //             },
    //         ).unwrap();
    //     })
    // });

    // c.bench_function("write_tiff_lu_bounded", |b| {
    //     b.iter(|| {
    //         write_tiff(
    //             &tile,
    //             &output_path,
    //             scale_shift,
    //             InterpolationStrategy::PlaneBased {
    //                 plane_interpolator,
    //                 plane_solver:lu_plane_solver,
    //             },
    //         ).unwrap();
    //     })
    // });

    // Clean up test output
    fs::remove_file(output_path).expect("Failed to clean up test output file");
}

criterion_group!(benches, benchmark);
criterion_main!(benches);