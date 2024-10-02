use qmlib::{quantized_mesh_tile, tiff_writer};
use std::path::PathBuf;
use std::env;
use std::result::Result; // Use standard Result

fn main() -> Result<(), String> {
    // Get the path from command line arguments
    let args: Vec<String> = env::args().collect();

    if args.len() < 3 {
        return Err(format!("Usage: {} <path> <scale shift>", args[0]));
    }
    let path_str = &args[1];
    let scale_shift = args[2].parse::<u16>().unwrap();
    let path: PathBuf = PathBuf::from(path_str);
    let mut outpath = path.clone();

    // Load the tile
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path)?;

    // write to the same location
    outpath.set_extension("tiff");
    tiff_writer::write_tiff(&tile, &outpath, scale_shift).map_err(|e| format!("Error exporting to KML: {}", e))?;


    Ok(())
}