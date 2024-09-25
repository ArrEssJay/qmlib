use core::str;
use qmlib::{tile, GeodeticPoint3, ToDegrees};
use std::path::PathBuf;
use std::env;
use std::result::Result; // Add this import to use Result

fn main() -> Result<(), String> { // Update Result to specify the error type
    // Get the path from command line arguments
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <path>", args[0]);
        std::process::exit(1);
    }

    let pathstr = &args[1];
    let path: PathBuf = PathBuf::from(pathstr);

    // Decode zoom, x, y from the file path
    match tile::load_quantized_mesh(&path) {
        Ok(tile) => {
            // Get bounding box
            println!("Tile Bounding Rectangle: {:#?}", tile.quantized_mesh.bounding_box);
            println!(
                "Tile Bounding Rectangle Dimensions: {:#?}",
                tile.quantized_mesh.bounding_box.dimensions()
            );

            // Process extensions - TODO
            for ext in &tile.quantized_mesh.extensions {
                println!(" - {:?}", ext.extension_id);
                if ext.extension_id == 4 {
                    let s = str::from_utf8(&ext.extension_data)
                        .unwrap_or_else(|e| panic!("Invalid UTF-8 sequence: {}", e));
                    println!("Metadata extension: \n{:?}", s);
                }
            }
            // Print center of tile in ecef and lat/lon
            println!("Tile Centre (Geocentric): {:#?}", tile.quantized_mesh.header.center);

            let centre_geodetic: GeodeticPoint3 =
                qmlib::ecef_to_geodetic(&tile.quantized_mesh.header.center, &tile.quantized_mesh.ellipsoid);
            println!(
                "Tile Centre (Geodetic): {:#?}",
                centre_geodetic.to_degrees()
            );

            // Calculate the ENU to ECEF rotation matrix and the distance to the bounding sphere center
            let to_enu_matrix =
                qmlib::calculate_enu_to_ecef_rotation_matrix(tile.quantized_mesh.header.center, &tile.quantized_mesh.ellipsoid)
                    .transpose();
            let dist_enu = to_enu_matrix
                .transform_vector(&(tile.quantized_mesh.header.bounding_sphere.center - tile.quantized_mesh.header.center));
            println!(
                "Tile Centre -> Bounding Sphere Centre (ENU from tile centre): {:#?}",
                dist_enu
            );

            // Print the first triangle
            println!("0: {:#?}", tile.quantized_mesh.get_triangle(0, true));

            // Print the first triangle
            println!("Vertices: {:#?}", tile.quantized_mesh.geodetic_vertices);

            println!("Vertices u: {:#?}", tile.quantized_mesh.vertex_data.u);
            println!("Vertices v: {:#?}", tile.quantized_mesh.vertex_data.v);
            println!("Vertices h: {:#?}", tile.quantized_mesh.vertex_data.height);
        }
        Err(err) => {
            eprintln!("Error: {}", err);
        }
    }

    Ok(())
}
