use std::sync::{atomic::{AtomicU32, Ordering}, Arc};

use nalgebra::{ Matrix2, Point2};
use rayon::prelude::*;

use crate::{geometry::Triangle, quantized_mesh_tile::QuantizedMeshTile, UV_SIZE_U16};

pub fn raster_dim_pixels(scale_shift: u16) -> u16 {
    UV_SIZE_U16 >> scale_shift
}

pub fn interpolate_triangle_face_height(
    point_usize: Point2<usize>,
    triangle_u16: &Triangle<u16>,
    heights: &[f64; 3],
) -> Option<f32> {

    // Cast triangle to f64 for bounds checking precision
    let triangle: Triangle<f64> = Triangle {
        vertices: triangle_u16.vertices.map(|vertex| {
            Point2::new(vertex.x as f64, vertex.y as f64)
        })
    };

    let point:Point2<f64> = Point2::new(point_usize.x as f64, point_usize.y as f64);  
    let v0 = triangle.vertices[0];
    let v1 = triangle.vertices[1];
    let v2 = triangle.vertices[2];

    // Translate points so that v0 is the origin
    let p_prime = point - v0;
    let v1_prime = v1 - v0;
    let v2_prime = v2 - v0;

    // Create the 2x2 matrix A with v1_prime and v2_prime as columns
    let a = Matrix2::new(v1_prime.x, v2_prime.x, v1_prime.y, v2_prime.y);

    // Check if the matrix is invertible (non-zero determinant)
    if let Some(a_inv) = a.try_inverse() {
        // Solve for lambda1 and lambda2
        let lambda = a_inv * p_prime;

        // Calculate lambda0
        let lambda0 = 1.0 - lambda.x - lambda.y;

        // Ensure that the barycentric coordinates are valid (inside the triangle)
        // Add a small epsilon threshold to account for floating-point precision errors.
        if lambda0 >= -f64::EPSILON && lambda.x >= -f64::EPSILON && lambda.y >= -f64::EPSILON {
            // Interpolate the height using the barycentric coordinates
            // Reduce precision to f32
            #[allow(clippy::cast_possible_truncation)]
            let interpolated_height = (lambda0 * heights[0]
                + lambda.x * heights[1] 
                + lambda.y * heights[2] ) as f32;

            Some(interpolated_height)
        } else {
            // The point is outside the triangle
            None
        }
    } else {
        // The triangle is degenerate (area is zero), no solution
        None
    }
}


pub fn rasterise(qmt: &QuantizedMeshTile, scale_shift: u16) -> Vec<f32> {
       
    let raster_dim_size:usize =  usize::from(raster_dim_pixels(scale_shift));
    let raster_size = raster_dim_size  * raster_dim_size ;
    let v = &qmt.quantized_mesh.vertex_data;
    let heights = &qmt.quantized_mesh.interpolated_height_vertices();

    // Store F32 in AtomicU32 using F32 bit patterns
    // Initialise with NAN so we can test whether a cell has already been written
    let raster: Vec<AtomicU32> = (0..raster_size as usize)
    .map(|_| AtomicU32::new(f32::NAN.to_bits()))
    .collect();

    // Thread shared
    let raster: Arc<Vec<AtomicU32>> = Arc::new(raster);

    // Use par_iter() to parallelize the processing of triangle indices
    v.triangle_index.par_iter().for_each(|triangle_vertices| {
        
        let vertex_at = |i: usize| Point2::new(
            v.u[triangle_vertices[i] as usize ] >> scale_shift,
            v.v[triangle_vertices[i] as usize ] >> scale_shift,
        );

        let triangle = Triangle {   
         vertices: [
                vertex_at(0),
                vertex_at(1),
                vertex_at(2),
            ]
        };

        let tri_heights = [
            heights[triangle_vertices[0] as usize],
            heights[triangle_vertices[1] as usize],
            heights[triangle_vertices[2] as usize],
        ];

        // 0,0 is in the lower left hand corner of the raster
        // scan the pixels in the raster column, reverse row order
        // now 0,0 is in the upper left hand corner of the raster
        // Why? GDAL does not honour the geotiff spec regarding 
        // raster <-> world axis convention flexibilty 

        let triangle_bounds = triangle.bounding_rect();
        
        // To reset the scanline x each row
        let raster_x_origin = triangle_bounds.lower_left.x as usize;
        
        // Invert y-axis by subtracting raster_point.y from raster_dim_size - 1
        let raster_y = (raster_dim_size - 1) - (triangle_bounds.lower_left.y as usize);

        // Build 2D point as its needed for vector math in the triangle rasteriser
        let mut raster_point:Point2<usize> = Point2::new(raster_x_origin,raster_y);
        
        // NaN bit pattern
        let f32_bits_nan = f32::to_bits(f32::NAN);

        for _ in 0..(triangle_bounds.upper_right.y -triangle_bounds.lower_left.y) // range is inclusive..exclusive
        {
            raster_point.x = raster_x_origin; // reset scanline x each row
            for _ in 0..(triangle_bounds.upper_right.x -triangle_bounds.lower_left.x)
            {   
                // index in the flat raster
                let raster_idx = (raster_y * raster_dim_size) + raster_point.x;
            
                // Check if the raster cell is empty by reading the atomic value without locking
                let cell = &raster[raster_idx];
                if cell.load(Ordering::Relaxed) == f32_bits_nan {
                    if let Some(interpolated_height) =
                        interpolate_triangle_face_height(raster_point, &triangle, &tri_heights)
                    {
                        // Try to write the value if the cell is still "empty" (contains NaN)
                        cell.store(interpolated_height.to_bits(), Ordering::Relaxed);
                    }
                }
                raster_point.x += 1;
            }
            raster_point.y -= 1;
        };
    });


    let image_data: Vec<f32> = raster
        .iter() // Use `iter()` to borrow the values
        .map(|cell| f32::from_bits(cell.load(Ordering::Relaxed)))
        .collect();
    
    #[cfg(debug_assertions)]
    {
        eprintln!("Warning: logging spurious rasteriser behaviour.");

        // Check for NaN values and issue a warning
        let mut nan_count = 0;    
        for (i, &pixel) in image_data.iter().enumerate() {
            if pixel.is_nan() {
                let x = i % raster_dim_size;
                let y = i / raster_dim_size;
                eprintln!("NaN found at position ({}, {})", x, y);
                nan_count += 1;
            }
        }

        if nan_count > 0 {
            eprintln!("Warning: {} pixels in the image data are NaN.", nan_count);
        }
    }

    image_data
}
