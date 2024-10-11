use crate::quantized_mesh_tile::QuantizedMeshTile;
use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

pub fn write_csv(qmt: &QuantizedMeshTile, file_path: &Path) -> io::Result<()> {

    let vertices = &qmt.vertices_as_geodetic_point3().iter().map(|v| v.to_ecef(&qmt.ellipsoid)).collect::<Vec<_>>();
       
    let mut csv_data = String::new();
    csv_data.push_str("x, y, z\n");

    for triangle in &qmt.quantized_mesh.vertex_data.triangle_index {
        for vi in triangle {
            let v = &vertices[*vi as usize];
            let x = (v.0.x * 100.0).round() / 100.0;
            let y = (v.0.y * 100.0).round() / 100.0;
            let z = (v.0.z * 100.0).round() / 100.0;
            csv_data.push_str(&format!("{:.2}, {:.2}, {:.2}\n", x, y, z));
        }
    }

    let mut file = File::create(file_path)?;
    

    file.write_all(csv_data.as_bytes())

}
