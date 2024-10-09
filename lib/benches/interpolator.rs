use criterion::{criterion_group, criterion_main, Criterion};
use nalgebra::Point2;
use qmlib::{geometry::Triangle, interpolator::{interpolate_height_barycentric, interpolate_height_lu, interpolate_height_parametric}};

fn benchmark(c: &mut Criterion) {
   
    //  Define a single triangle
    let triangle_u16 = Triangle {
        vertices: [
            Point2::new(0, 0),
            Point2::new(10, 0),
            Point2::new(0, 10),
        ],
        };
        let heights = [0.0, 10.0, 20.0];
    
        // Create a grid of points that overlap the triangle
        let points: Vec<Point2<usize>> = (0..11)
            .flat_map(|x| (0..11).map(move |y| Point2::new(x, y)))
            .collect();
        
    
        c.bench_function("interpolate_height_barycentric", |b| {
            b.iter(|| {
                for point in &points {
                    let _ = interpolate_height_barycentric(*point, &triangle_u16, &heights);
                }
            })
        });
    
        c.bench_function("interpolate_height_lu", |b| {
            b.iter(|| {
                for point in &points {
                    let _ = interpolate_height_lu(*point, &triangle_u16, &heights, true); // Set bounds_check to true
                }
            })
        });
    
        c.bench_function("interpolate_height_parametric", |b| {
            b.iter(|| {
                for point in &points {
                    let _ = interpolate_height_parametric(*point, &triangle_u16, &heights);
                }
            })
        });

    // Clean up test output
}

criterion_group!(benches, benchmark);
criterion_main!(benches);