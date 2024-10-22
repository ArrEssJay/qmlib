use qmlib::interpolator::InterpolationMethod;
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

    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path)?;

    // write to the same location
    outpath.set_extension("tiff");
    tiff_writer::write_tiff(&tile, &outpath, scale_shift, InterpolationMethod::Barycentric).map_err(|e| format!("Error exporting to GeoTIFF: {}", e))?;


    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::fs;

    #[test]
    fn test_write_tiff() {
        // Assume you have a valid test tile path
        let mut path: PathBuf = PathBuf::from(format!("{}/../test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));       
        let scale_shift: u16 = 4; 
        let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();
        path.set_extension(".test.tiff");
        let result = tiff_writer::write_tiff(&tile, &path, scale_shift, InterpolationMethod::Barycentric);
        assert!(result.is_ok(), "Failed to write TIFF: {:?}", result.err());
        
        // Clean up test output
        fs::remove_file(path).expect("Failed to clean up test output file");
    }
}
