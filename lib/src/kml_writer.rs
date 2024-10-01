use crate::geometry::GeodeticPoint3;
use crate::QuantizedMesh;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn export_to_kml(qm: &QuantizedMesh, file_path: &Path) -> io::Result<()> {
    let mut kml_data = String::new();

    // KML Header
    kml_data.push_str(r#"<?xml version="1.0" encoding="UTF-8"?>"#);
    kml_data.push_str(r#"<kml xmlns="http://www.opengis.net/kml/2.2">"#);
    kml_data.push_str(r#"<Document>"#);

    // Get index data
    let index_data: &Vec<[usize; 3]> = &qm.vertex_data.triangle_index;
    let vertices = &qm.vertices_as_geodetic_point3();

    // Loop through triangles and write to KML
    for (idx, triangle_indices) in index_data.iter().enumerate() {
        
        let triangle: [&GeodeticPoint3<f64>; 3] = [
            &vertices[triangle_indices[0]],
            &vertices[triangle_indices[1]],
            &vertices[triangle_indices[2]],
        ];

        // KML Placemark for each triangle
        kml_data.push_str(r#"<Placemark>"#);
        kml_data.push_str(r#"<name>Triangle "#);
        kml_data.push_str(&idx.to_string());
        kml_data.push_str(r#"</name>"#);
        kml_data.push_str(r#"<Polygon><outerBoundaryIs><LinearRing><coordinates>"#);

        for point in triangle.iter() {
            // Convert each point to KML coordinates (lon, lat, alt)
            kml_data.push_str(&format!(
                "{},{},{}\n",
                point.lon(),
                point.lat(),
                point.alt()
            ));
        }

        kml_data.push_str(r#"</coordinates></LinearRing></outerBoundaryIs></Polygon></Placemark>"#);
    }

    // Add all geodetic vertices as placemarks
    for (i, vertex) in vertices.iter().enumerate() {
        kml_data.push_str(r#"<Placemark>"#);
        kml_data.push_str(&format!(r#"<name>Vertex {}</name>"#, i));
        kml_data.push_str(r#"<Point><coordinates>"#);
        kml_data.push_str(&format!(
            "{},{},{}\n",
            vertex.lon(),
            vertex.lat(),
            vertex.alt()
        ));
        kml_data.push_str(r#"</coordinates></Point></Placemark>"#);
    }

    // KML Footer
    kml_data.push_str(r#"</Document></kml>"#);

    // Write to file
    let mut file = File::create(file_path)?;
    file.write_all(kml_data.as_bytes())?;

    Ok(())
}
