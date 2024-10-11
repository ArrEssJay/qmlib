use core::str;
use geo::{Coord, GeodesicArea, LineString, Polygon};
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
                "Tile Bounding Rectangle: lat:{:#?}° lon:{:#?}°",
                qmt.bounding_rectangle.lower_left.lat().to_degrees(),
                qmt.bounding_rectangle.lower_left.lon().to_degrees()
            );

            let centre_geodetic = qmt.quantized_mesh.header.center.to_geodetic(&qmt.ellipsoid);
            println!(
                "Tile Centre (Geodetic): {:#?}",
                centre_geodetic.to_degrees()
            );
            // Print center of tile in ecef and lat/lon
            println!(
                "Tile Centre (Geocentric): {:#?}",
                qmt.quantized_mesh.header.center
            );

            // Print center of tile in ecef and lat/lon
            println!(
                "Tile Bounding Sphere Centre (Geocentric): {:#?}",
                qmt.quantized_mesh.header.bounding_sphere.center
            );
            println!(
                "Tile Bounding Sphere Radius: {:#?}m",
                qmt.quantized_mesh.header.bounding_sphere.radius
            );

            // Determine tile centre (ellipsoid height)
            // and bounding sphere centre offset
            // For validation these should be a distance vertically
            // approximately equal to the average tile height
            //
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

            println!(
                "Tile Height min: {:#?}, max: {:#?}",
                &qmt.quantized_mesh.header.min_height,
                &qmt.quantized_mesh.header.max_height
            );

            let br = qmt.bounding_rectangle;
            let min_x = br.lower_left.0.x;
            let max_x = br.upper_right.0.x;
            let min_y = br.lower_left.0.y;
            let max_y = br.upper_right.0.y;

            let exterior = LineString::from(vec![
                Coord { x: min_x, y: min_y },
                Coord { x: min_x, y: max_y },
                Coord { x: max_x, y: max_y },
                Coord { x: max_x, y: min_y },
            ]);

            // Determine the tile area
            let polygon = Polygon::new(exterior, vec![]);

            let area = polygon.geodesic_area_unsigned();
            let formatted_area = format!("{:.2e}", area);
            println!("\nBounding Rectangle Area: {formatted_area}m²");

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
