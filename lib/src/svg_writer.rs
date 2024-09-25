// svg_writer.rs
use crate::QuantizedMesh;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn export_to_svg(quantized_mesh: &QuantizedMesh, file_path: &Path) -> io::Result<()> {
    let mut svg_data = String::new();

    // SVG Header
    svg_data.push_str(r#"<?xml version="1.0" standalone="no"?>"#);
    svg_data.push_str(
        r#"<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="32767" height="32767">"#,
    );

    let index_data = quantized_mesh
        .vertex_data
        .index_data()
        .expect("Failed to retrieve index data");

        for triangle_index in 0..(index_data.len() / 3) {
            let tri = quantized_mesh.get_triangle(triangle_index, false)
                .expect(&format!("triangle at index {} could not be retrieved", triangle_index));
        
            // Create SVG polygon element
            svg_data.push_str(r#"<polygon points=""#);
            
            // Collect points as a single string
            let points: String = tri.iter()
                .map(|point| format!("{} {}", point.x, point.y)) // Assuming Point3<f64> with x and y
                .collect::<Vec<String>>()
                .join(" "); // Join points with a space
        
            // Append the points directly, without additional quotes around them
            svg_data.push_str(&points); // Append points
            svg_data.push_str(r#"" style="fill:lime;stroke:purple;stroke-width:1" />"#); // Close polygon element
        }

    // SVG Footer
    svg_data.push_str(r#"</svg>"#);

    // Write to file
    let mut file = File::create(file_path)?;
    file.write_all(svg_data.as_bytes())?;

    Ok(())
}
