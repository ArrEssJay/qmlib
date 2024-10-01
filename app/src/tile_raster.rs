use qmlib::{tile, tiff_writer};
use std::path::PathBuf;
use std::env;
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
    let tile = tile::load_quantized_mesh(&path)?;

    let mut outpath = path.clone();
        
    outpath.set_extension("tiff");
    tiff_writer::export_to_tiff(&tile.quantized_mesh, &outpath, 4).map_err(|e| format!("Error exporting to KML: {}", e))?;


    Ok(())
}