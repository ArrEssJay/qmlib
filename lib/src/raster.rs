use std::sync::{atomic::{AtomicU32, Ordering}, Arc};

use nalgebra:: Point2;
use rayon::prelude::*;

use crate::{geometry::Triangle, interpolator::Interpolator, quantized_mesh_tile::QuantizedMeshTile, UV_SIZE_U16};



pub fn raster_dim_pixels(scale_shift: u16) -> u16 {
    UV_SIZE_U16 >> scale_shift
}

pub fn rasterise(qmt: &QuantizedMeshTile, scale_shift: u16,  interpolator: Interpolator) -> Vec<f32> {
       
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
                        interpolator(raster_point, &triangle, &tri_heights)
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
