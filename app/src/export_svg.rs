use qmlib::{svg_writer, quantized_mesh_tile};
use std::env;
use std::path::PathBuf;
use std::result::Result; // Use standard Result

fn main() -> Result<(), String> {
    // Get the path from command line arguments
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        return Err(format!("Usage: {} <path>", args[0]));
    }

    let pathstr = &args[1];
    let path: PathBuf = PathBuf::from(pathstr);

    // Load the QuantizedMesh
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path)?;

    // Export to KML
    let mut outpath = path.clone();

    outpath.set_extension("svg");
    svg_writer::write_svg(&tile.quantized_mesh, &outpath)
        .map_err(|e| format!("Error exporting to SVG: {}", e))?;

    Ok(())
}
