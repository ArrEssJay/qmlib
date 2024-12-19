use clap::{Arg, Command};
use qmlib::{quantized_mesh_tile, tiff_writer};
use std::path::PathBuf;
use std::result::Result; // Use standard Result

fn main() -> Result<(), String> {
    let matches = Command::new("export_geotiff")
        .about("Exports a Quantized Mesh Tile to GeoTIFF")
        .arg(
            Arg::new("input")
                .help("Input path to the Quantized Mesh Tile")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::new("raster_scale_factor")
                .help("Scale shift value as an integer")
                .required(true)
                .index(2),
        )
        .arg(
            Arg::new("rasteriser")
                .help("Rasteriser to use: 'cpu' or 'gpu' (defaults to 'gpu')")
                .required(false)
                .index(3),
        )
        .get_matches();

    // Retrieve the input path
    let input_path: PathBuf = matches
        .get_one::<String>("input")
        .expect("Input path is required")
        .into();

    // Parse the raster_scale_factor argument
    let raster_scale_factor: u32 = matches
        .get_one::<String>("raster_scale_factor")
        .expect("Scale shift is required")
        .parse()
        .map_err(|_| "Invalid raster_scale_factor value")?;

    // Determine the rasteriser to use
    let rasteriser = match matches
        .get_one::<String>("rasteriser")
        .unwrap_or(&"gpu".to_string())
        .to_lowercase()
        .as_str()
    {
        "cpu" => tiff_writer::Rasteriser::CPU,
        "gpu" => tiff_writer::Rasteriser::GPU,
        _ => {
            eprintln!("Invalid rasteriser. Expected 'cpu' or 'gpu'.");
            return Err("Invalid rasteriser".to_string());
        }
    };

    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&input_path)?;

    let mut outpath = input_path.clone();
    outpath.set_extension("tiff");

    tiff_writer::write_tiff(&tile, &outpath, raster_scale_factor, rasteriser)
        .map_err(|e| format!("Error exporting to GeoTIFF: {}", e))?;

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;
    use std::fs;

    #[test]
    fn test_write_tiff() {
        // Assume you have a valid test tile path
        let mut path: PathBuf = PathBuf::from(format!(
            "{}/../test/terrain_data/a/15/59489/9692.terrain",
            env!("CARGO_MANIFEST_DIR")
        ));
        let raster_scale_factor: u32 = 4;
        let tile = quantized_mesh_tile::load_quantized_mesh_tile(&path).unwrap();
        path.set_extension(".test.tiff");
        let result =
            tiff_writer::write_tiff(&tile, &path, raster_scale_factor, tiff_writer::Rasteriser::CPU);
        assert!(result.is_ok(), "Failed to write TIFF: {:?}", result.err());

        // Clean up test output
        fs::remove_file(path).expect("Failed to clean up test output file");
    }
}
