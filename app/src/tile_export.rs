use binrw::{
    io::{Error, ErrorKind, Result},
    BinRead,
};
use qmlib::{Ellipsoid, QuantizedMesh};
use std::fs::File;
use std::path::PathBuf;


use qmlib::kml_writer;
use qmlib::svg_writer;

fn main() -> Result<()> {
    // Parameters for the file path
    let zoom = 18;
    let x = 470053;
    let y = 77832;
    let bounding_box = qmlib::tile_to_bbox_crs84(x, y, zoom);
    let ellipsoid: qmlib::Ellipsoid = Ellipsoid::wgs84();
    let tiling_scheme = qmlib::TilingScheme::Slippy;

    let pathstr = "./test/data/a";
    // Construct the file path
    let mut path: PathBuf = PathBuf::from(pathstr);
    path.push(zoom.to_string());
    path.push(x.to_string());
    path.push(format!("{}.terrain", y));
    println!("File: {:?}", path);

    // Open the file
    let mut file = File::open(&path)?;

    // Read the QuantizedMesh
    let qmdata: QuantizedMesh = QuantizedMesh::read_le(&mut file).map_err(|e| {
        Error::new(
            ErrorKind::Other,
            format!("Failed to read quantized mesh: {:?}", e),
        )
    })?;

    let qm = QuantizedMesh::new(
        qmdata.header,
        qmdata.vertex_data,
        qmdata.extensions,
        bounding_box,
        ellipsoid,
        tiling_scheme,
    );

    let kml_filename = format!("{}.kml", y);
    let mut outpath: PathBuf = PathBuf::from("./test/data/a");

    outpath.push(&kml_filename);
    if let Err(e) = kml_writer::export_to_kml(&qm, &outpath) {
        eprintln!("Error exporting to KML: {}", e);
    } else {
        println!("KML exported successfully to {}", &kml_filename);
    }


    let svg_filename = format!("{}.svg", y);
    let mut outpath: PathBuf = PathBuf::from("./test/data/a");

    outpath.push(&svg_filename);
    if let Err(e) = svg_writer::export_to_svg(&qm, &outpath) {
        eprintln!("Error exporting to SVG: {}", e);
    } else {
        println!("SVG exported successfully to {}", &svg_filename);
    }

    Ok(())
}
