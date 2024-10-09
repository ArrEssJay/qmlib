use std::env;
use nalgebra::Point2;
use qmlib::{geometry::GeodeticPoint2, quantized_mesh_tile::{tiles_for_point, CRS}};

fn main() {
    // Collect command-line arguments
    let args: Vec<String> = env::args().collect();

    // Check for correct number of arguments
    if args.len() < 3 {
        eprintln!("Usage: {} <lat> <lon> <zoom> [<CRS>]", args[0]);
        return;
    }

    // Parse latitude and longitude from arguments
    let lat: f64 = args[1].parse().expect("Invalid latitude");
    let lon: f64 = args[2].parse().expect("Invalid longitude");

    // Parse zoom level or default to 20
    let zoom: u8 = if args.len() > 3 {
        args[3].parse().expect("Invalid zoom level")
    } else {
        20
    };

    // Parse CRS from command-line arguments, defaulting to Epsg4326
    let crs = if args.len() > 4 {
        match args[4].to_lowercase().as_str() {
            "epsg3857" => CRS::Epsg3857,
            _ => CRS::Epsg4326,
        }
    } else {
        CRS::Epsg4326
    };

    // Create the GeodeticPoint2 instance
    let point = GeodeticPoint2(Point2::new(lon, lat)); // Note the order (lon, lat)

    // Call the tiles_for_point function
    match tiles_for_point(point, zoom, &crs) {
        Ok(tiles) => {
            // Print the results in a table format
            println!("{:<8} {:<8} {:<5}", "X", "Y", "Zoom");
            for (x, y, z) in tiles {
                println!("{:<8} {:<8} {:<5}", x, y, z);
            }
        }
        Err(e) => eprintln!("Error: {}", e),
    }
}
