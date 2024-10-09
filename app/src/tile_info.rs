use core::str;
use std::env;
use std::path::PathBuf;
use std::result::Result;

use qmlib::geometry::calculate_enu_to_ecef_rotation_matrix;
use qmlib::quantized_mesh_tile;

fn main() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        eprintln!("Usage: {} <path>", args[0]);
        std::process::exit(1);
    }

    let pathstr = &args[1];
    let path: PathBuf = PathBuf::from(pathstr);

    // Decode zoom, x, y from the file path
    match quantized_mesh_tile::load_quantized_mesh_tile(&path) {
        Ok(qmt) => {
            // Get bounding box
            println!(
                "Tile Bounding Rectangle: {:#?}",
                qmt.bounding_rectangle
            );

            // Print center of tile in ecef and lat/lon
            println!(
                "Tile Centre (Geocentric): {:#?}",
                qmt.quantized_mesh.header.center
            );

            let centre_geodetic = qmt
                .quantized_mesh
                .header
                .center
                .to_geodetic(&qmt.ellipsoid);
            println!(
                "Tile Centre (Geodetic): {:#?}",
                centre_geodetic.to_degrees()
            );

            //Calculate the ENU to ECEF rotation matrix and the distance to the bounding sphere center
            let to_enu_matrix = calculate_enu_to_ecef_rotation_matrix(
                &qmt.quantized_mesh.header.center,
                &qmt.ellipsoid,
            )
            .transpose();
            let dist_enu = to_enu_matrix.transform_vector(
                &(qmt.quantized_mesh.header.bounding_sphere.center.0
                    - qmt.quantized_mesh.header.center.0),
            );
            println!(
                "Tile Centre -> Bounding Sphere Centre (ENU from tile centre): {:#?}",
                dist_enu
            );

            println!("\nVertex Array: count => [min:max] ");
            println!(
                "Vertex (U): {:#?} => [{:#?}:{:#?}]",
                qmt.quantized_mesh.vertex_data.u.len(),
                qmt.quantized_mesh.vertex_data.u.iter().min().unwrap(),
                qmt.quantized_mesh.vertex_data.u.iter().max().unwrap()
            );
            println!(
                "Vertex (V): {:#?} => [{:#?}:{:#?}]",
                qmt.quantized_mesh.vertex_data.v.len(),
                qmt.quantized_mesh.vertex_data.v.iter().min().unwrap(),
                qmt.quantized_mesh.vertex_data.v.iter().max().unwrap()
            );
            println!(
                "Vertex (H): {:#?} => [{:#?}:{:#?}]",
                qmt.quantized_mesh.vertex_data.height.len(),
                qmt.quantized_mesh.vertex_data.height.iter().min().unwrap(),
                qmt.quantized_mesh.vertex_data.height.iter().max().unwrap()
            );

            println!(
                "\n\nTriangle Count: {:#?}",
                qmt.quantized_mesh.vertex_data.triangle_count
            );
            println!(
                "East Edge Vertex Count: {:?}",
                qmt.quantized_mesh
                    .vertex_data
                    .edge_indices
                    .east_vertex_count
            );
            println!(
                "West Edge Vertex Count: {:?}",
                qmt.quantized_mesh
                    .vertex_data
                    .edge_indices
                    .west_vertex_count
            );
            println!(
                "North Edge Vertex Count: {:?}",
                qmt.quantized_mesh
                    .vertex_data
                    .edge_indices
                    .north_vertex_count
            );
            println!(
                "South Edge Vertex Count: {:?}",
                qmt.quantized_mesh
                    .vertex_data
                    .edge_indices
                    .south_vertex_count
            );

            // Process extensions - TODO
            for ext in &qmt.quantized_mesh.extensions {
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
