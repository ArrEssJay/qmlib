// svg_writer.rs
use crate::{QuantizedMesh, VertexData};
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn write_svg(quantized_mesh: &QuantizedMesh, file_path: &Path) -> io::Result<()> {
    let mut svg_data = String::new();

    // SVG Header
    svg_data.push_str(r#"<?xml version="1.0" standalone="no"?>"#);
    svg_data.push_str(
        r#"<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="32767" height="32767">"#,
    );

    let vertex_data: &VertexData = &quantized_mesh.vertex_data;
    let triangle_index = &vertex_data.triangle_index;

    for vertex in triangle_index {
        let u0 = vertex_data.u[vertex[0] as usize];
        let v0 = vertex_data.v[vertex[0] as usize];

        let u1 = vertex_data.u[vertex[1] as usize];
        let v1 = vertex_data.v[vertex[1] as usize];

        let u2 = vertex_data.u[vertex[2] as usize];
        let v2 = vertex_data.v[vertex[2] as usize];
        // Create SVG polygon element
        svg_data.push_str(r#"<polygon points=""#);

        // Collect points as a single string
        let points: String = [
            format!("{u0} {v0}"),
            format!("{u1} {v1}"),
            format!("{u2} {v2}"),
        ]
        .join(" "); // Join points with a space

        // Append the points directly, without additional quotes around them
        svg_data.push_str(&points); // Append points
        svg_data.push_str(r#"" style="fill:lime;stroke:purple;stroke-width:1" />"#);
        // Close polygon element
    }

    // SVG Footer
    svg_data.push_str(r"</svg>");

    // Write to file
    let mut file = File::create(file_path)?;
    file.write_all(svg_data.as_bytes())?;

    Ok(())
}
