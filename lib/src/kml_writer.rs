use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use crate::QuantizedMesh;


    pub fn export_to_kml(quantized_mesh: &QuantizedMesh, file_path: &Path) -> io::Result<()> {
        let mut kml_data = String::new();

        // KML Header
        kml_data.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
        kml_data.push_str(r#"<kml xmlns="http://www.opengis.net/kml/2.2">"#);
        kml_data.push_str(r#"<Document>"#);

        // Get index data
        let index_data = quantized_mesh
            .vertex_data
            .index_data()
            .expect("Failed to retrieve index data");

        // Loop through triangles and write to KML
        for triangle_index in 0..(index_data.len() / 3) {
            // Use expect to unwrap the Option
            let triangle = quantized_mesh.get_triangle(triangle_index, true).expect(&format!(
                "Triangle at index {} could not be retrieved",
                triangle_index
            ));

            // KML Placemark for each triangle
            kml_data.push_str(r#"<Placemark>"#);
            kml_data.push_str(r#"<name>Triangle "#);
            kml_data.push_str(&triangle_index.to_string());
            kml_data.push_str(r#"</name>"#);
            kml_data.push_str(r#"<Polygon><outerBoundaryIs><LinearRing><coordinates>"#);

            for point in triangle.iter() {
                // Convert each point to KML coordinates (lon, lat, alt)
                kml_data.push_str(&format!("{},{},{}\n", point.x, point.y, point.z));
            }

            kml_data
                .push_str(r#"</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>"#);
        }

        // Add all geodetic vertices as placemarks
        for (i, vertex) in quantized_mesh.geodetic_vertices.iter().enumerate() {
            kml_data.push_str(r#"<Placemark>"#);
            kml_data.push_str(&format!(r#"<name>Vertex {}</name>"#, i));
            kml_data.push_str(r#"<Point><coordinates>"#);
            kml_data.push_str(&format!("{},{},{}\n", vertex[0], vertex[1], vertex[2]));
            kml_data.push_str(r#"</coordinates></Point></Placemark>"#);
        }


        // KML Footer
        kml_data.push_str(r#"</Document></kml>"#);

        // Write to file
        let mut file = File::create(file_path)?;
        file.write_all(kml_data.as_bytes())?;

        Ok(())
    }
