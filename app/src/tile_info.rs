use core::str;
use qmlib::{tile, GeodeticPoint3, ToDegrees};
use std::env;
use std::path::PathBuf;
use std::result::Result; // Add this import to use Result

fn main() -> Result<(), String> {
    // Update Result to specify the error type
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
            println!(
                "Tile Bounding Rectangle: {:#?}",
                tile.quantized_mesh.bounding_box
            );

            // Print center of tile in ecef and lat/lon
            println!(
                "Tile Centre (Geocentric): {:#?}",
                tile.quantized_mesh.header.center
            );

            let centre_geodetic: GeodeticPoint3 = qmlib::ecef_to_geodetic(
                &tile.quantized_mesh.header.center,
                &tile.quantized_mesh.ellipsoid,
            );
            println!(
                "Tile Centre (Geodetic): {:#?}",
                centre_geodetic.to_degrees()
            );

            // Calculate the ENU to ECEF rotation matrix and the distance to the bounding sphere center
            let to_enu_matrix = qmlib::calculate_enu_to_ecef_rotation_matrix(
                tile.quantized_mesh.header.center,
                &tile.quantized_mesh.ellipsoid,
            )
            .transpose();
            let dist_enu = to_enu_matrix.transform_vector(
                &(tile.quantized_mesh.header.bounding_sphere.center
                    - tile.quantized_mesh.header.center),
            );
            println!(
                "Tile Centre -> Bounding Sphere Centre (ENU from tile centre): {:#?}",
                dist_enu
            );

            println!(
                "\n\nVertex Count (U): {:#?}",
                tile.quantized_mesh.vertex_data.u.len()
            );
            println!(
                "Vertex Count (V): {:#?}",
                tile.quantized_mesh.vertex_data.v.len()
            );
            println!(
                "Vertex Count (Height): {:#?}",
                tile.quantized_mesh.vertex_data.height.len()
            );

            println!(
                "\n\nTriangle Count: {:#?}",
                tile.quantized_mesh.vertex_data.triangle_count
            );
            println!(
                "Triangle Long Index: {:#?}",
                tile.quantized_mesh.vertex_data.index_data_long.is_some()
            );
            println!(
                "Triangle Short Index: {:#?}",
                tile.quantized_mesh.vertex_data.index_data_short.is_some()
            );

            println!(
                "\nLong Edge Index: {:#?}",
                tile.quantized_mesh.vertex_data.edge_indices_long.is_some()
            );
            println!(
                "Short Edge Index: {:#?}",
                tile.quantized_mesh.vertex_data.edge_indices_short.is_some()
            );

            if let Some(edge_indices_long) = &tile.quantized_mesh.vertex_data.edge_indices_long {
                // Now you can use edge_indices_long safely
                // For example:
                println!(
                    "East Edge Vertex Count: {:?}",
                    edge_indices_long.east_vertex_count
                );
                println!(
                    "West Edge Vertex Count: {:?}",
                    edge_indices_long.west_vertex_count
                );
                println!(
                    "North Edge Vertex Count: {:?}",
                    edge_indices_long.north_vertex_count
                );
                println!(
                    "South Edge Vertex Count: {:?}",
                    edge_indices_long.south_vertex_count
                );
            }

            if let Some(edge_indices_short) = &tile.quantized_mesh.vertex_data.edge_indices_short {
                // Now you can use edge_indices_long safely
                // For example:
                println!(
                    "East Edge Vertex Count: {:?}",
                    edge_indices_short.east_vertex_count
                );
                println!(
                    "West Edge Vertex Count: {:?}",
                    edge_indices_short.west_vertex_count
                );
                println!(
                    "North Edge Vertex Count: {:?}",
                    edge_indices_short.north_vertex_count
                );
                println!(
                    "South Edge Vertex Count: {:?}",
                    edge_indices_short.south_vertex_count
                );
            }

            // Process extensions - TODO
            for ext in &tile.quantized_mesh.extensions {
                if ext.extension_id == 4 {
                    let s = str::from_utf8(&ext.extension_data)
                        .unwrap_or_else(|e| panic!("Invalid UTF-8 sequence: {}", e));
                    println!("\nMetadata Extension (4): \n{:?}", s);
                }
            }
        }
        Err(err) => {
            eprintln!("Error: {}", err);
        }
    }

    Ok(())
}
