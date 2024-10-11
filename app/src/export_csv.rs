use qmlib::{csv_writer, quantized_mesh_tile};
use std::env;
use std::path::PathBuf;
use std::result::Result; // Use standard Result

fn main() -> Result<(), String> {
    // Get the path from command line arguments
    let args: Vec<String> = env::args().collect();

    if args.len() < 2 {
        return Err(format!("Usage: {} <path> [--geodetic|--geocentric]", args[0]));
    }

    let mode = if args.len() > 2 {
        match args[1].as_str() {
            "--geodetic" => "geodetic",
            "--geocentric" => "geocentric",
            _ => return Err(format!("Invalid option: {}. Use --geodetic (default) or --geocentric.", args[1])),
        }
    } else {
        "geodetic" // default mode
    };

    let coord_type = if mode == "geodetic" {
        csv_writer::CoordType::Geodetic
    } else {
        csv_writer::CoordType::Geocentric
    };

    println!("Outputting {} co-ordinates", mode);

    let pathstr = &args[2];
    let path: PathBuf = PathBuf::from(pathstr);

    // Load the QuantizedMesh
    let qmt = quantized_mesh_tile::load_quantized_mesh_tile(&path)?;

    // Export to KML
    let mut outpath = path.clone();

    outpath.set_extension("csv");
    csv_writer::write_csv(&qmt, &outpath, coord_type)
        .map_err(|e| format!("Error exporting to CSV: {}", e))?;

    Ok(())
}
