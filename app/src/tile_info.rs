use binrw::{
    io::{Error, ErrorKind, Result},
    BinRead,
};
use core::str;
use qmlib::{Ellipsoid, GeodeticPoint3, QuantizedMesh, ToDegrees};
use std::fs::File;
use std::path::PathBuf;

fn main() -> Result<()> {
    // Parameters for the file path
    let zoom = 18;
    let x = 470053;
    let y = 77832;
    let bounding_box = qmlib::tile_to_bbox_crs84(x, y, zoom);
    let ellipsoid: qmlib::Ellipsoid = Ellipsoid::wgs84();
    let tiling_scheme = qmlib::TilingScheme::Slippy;

    // Construct the file path
    let mut path = PathBuf::from("./test/data/a");
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

    // Get bounding box
    println!("Tile Bounding Rectangle: {:#?}", qm.bounding_box);
    println!(
        "Tile Bounding Rectangle Dimensions: {:#?}",
        qm.bounding_box.dimensions()
    );

    // Process extensions - TODO
    for ext in &qm.extensions {
        println!(" - {:?}", ext.extension_id);
        if ext.extension_id == 4 {
            let s = str::from_utf8(&ext.extension_data)
                .unwrap_or_else(|e| panic!("Invalid UTF-8 sequence: {}", e));
            println!("Metadata extension: \n{:?}", s);
        }
    }
    // Print center of tile in ecef and lat/lon
    println!("Tile Centre (Geocentric): {:#?}", qm.header.center);

    let centre_geodetic: GeodeticPoint3 = qmlib::ecef_to_geodetic(&qm.header.center, &qm.ellipsoid);
    println!(
        "Tile Centre (Geodetic): {:#?}",
        centre_geodetic.to_degrees()
    );

    // Calculate the ENU to ECEF rotation matrix and the distance to the bounding sphere center
    let to_enu_matrix =
        qmlib::calculate_enu_to_ecef_rotation_matrix(qm.header.center, &qm.ellipsoid).transpose();
    let dist_enu =
        to_enu_matrix.transform_vector(&(qm.header.bounding_sphere.center - qm.header.center));
    println!(
        "Tile Centre -> Bounding Sphere Centre (ENU from tile centre): {:#?}",
        dist_enu
    );

    // Print the first triangle
    println!("0: {:#?}", qm.get_triangle(0, &qm.bounding_box));

    Ok(())
}
