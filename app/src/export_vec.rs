use qmlib::{kml_writer, svg_writer, twodm_writer, quantized_mesh_tile};
use std::env;
use std::path::PathBuf;
use std::result::Result; // Use standard Result
use std::process;

fn main() -> Result<(), String> {
    // Get the arguments from command line
    let args: Vec<String> = env::args().collect();

    // Check for correct usage
    if args.len() < 3 {
        eprintln!("Usage: {} <input_path> <output_format>", args[0]);
        eprintln!("Supported output formats: kml, svg, 2dm");
        return Err("Invalid arguments".to_string());
    }

    let input_path_str = &args[1];
    let output_format = &args[2].to_lowercase(); // Normalize to lowercase
    let input_path: PathBuf = PathBuf::from(input_path_str);

    // Load the QuantizedMeshTile
    let tile = quantized_mesh_tile::load_quantized_mesh_tile(&input_path)
        .map_err(|e| format!("Error loading quantized mesh tile: {}", e))?;

    // Determine the output path and writer function based on the format
    let (output_path, writer_fn): (PathBuf, fn(&quantized_mesh_tile::QuantizedMeshTile, &PathBuf) -> Result<(), String>) = match output_format.as_str() {
        "kml" => {
            let mut outpath = input_path.clone();
            outpath.set_extension("kml");
            (outpath, kml_write_wrapper)
        },
        "svg" => {
            let mut outpath = input_path.clone();
            outpath.set_extension("svg");
            (outpath, svg_write_wrapper)
        },
        "2dm" => {
            let mut outpath = input_path.clone();
            outpath.set_extension("2dm");
            (outpath, twodm_write_wrapper)
        },
        _ => {
            eprintln!("Unsupported output format: {}", output_format);
            eprintln!("Supported formats are: kml, svg, 2dm");
            process::exit(1);
        }
    };

    // Call the appropriate writer function
    writer_fn(&tile, &output_path)?;

    println!("Exported to {}", output_path.display());

    Ok(())
}

// Wrapper functions to match the expected function type
fn kml_write_wrapper(tile: &quantized_mesh_tile::QuantizedMeshTile, path: &PathBuf) -> Result<(), String> {
    kml_writer::write_kml(tile, path)
        .map_err(|e| format!("Error exporting to KML: {}", e))
}

fn svg_write_wrapper(tile: &quantized_mesh_tile::QuantizedMeshTile, path: &PathBuf) -> Result<(), String> {
    svg_writer::write_svg(&tile.quantized_mesh, path)
        .map_err(|e| format!("Error exporting to SVG: {}", e))
}

fn twodm_write_wrapper(tile: &quantized_mesh_tile::QuantizedMeshTile, path: &PathBuf) -> Result<(), String> {
    twodm_writer::write_2dm(tile, path)
        .map_err(|e| format!("Error exporting to 2DM: {}", e))
}