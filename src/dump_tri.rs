use binrw::{
    io::{Error, ErrorKind, Result},
    BinRead,
};
use core::str;
use std::fs::File;
use std::path::PathBuf;
mod quantized_mesh;

fn main() -> Result<()> {
    // Parameters for the file path
    let zoom = 18;
    let x = 470053;
    let y = 77832;

    // Construct the file path
    let mut path = PathBuf::from("./test/data/a");
    path.push(zoom.to_string());
    path.push(x.to_string());
    path.push(format!("{}.terrain", y));
    println!("File: {:?}", path);

    // Open the file
    let mut file = File::open(&path)?;

    // Read the QuantizedMesh
    let qm: quantized_mesh::QuantizedMesh = quantized_mesh::QuantizedMesh::read(&mut file)
        .map_err(|e| {
            Error::new(
                ErrorKind::Other,
                format!("Failed to read quantized mesh: {:?}", e),
            )
        })?;

    // Process extensions
    for ext in &qm.extensions {
        println!(" - {:?}", ext.extensionId);
        if ext.extensionId == 4 {
            let s = str::from_utf8(&ext.extensionData)
                .unwrap_or_else(|e| panic!("Invalid UTF-8 sequence: {}", e));
            println!("Metadata extension: \n{:?}", s);
        }
    }

    // Get bounding box
    let bbox = quantized_mesh::tile_to_bbox_crs84(x, y, zoom);
    println!("Tile Bounding Rectangle: {:#?}", bbox);
    println!(
        "Tile Bounding Rectangle Dimensions: {:#?}",
        bbox.dimensions()
    );

    // Print center of tile in lat/lon
    let centre_geodetic =
        quantized_mesh::ecef_to_geodetic(&qm.header.center, &quantized_mesh::Ellipsoid::wgs84());
    println!("Tile Centre (Geodetic): {:#?}", centre_geodetic);

    // Calculate the ENU to ECEF rotation matrix and the distance to the bounding sphere center
    let to_enu_matrix = quantized_mesh::calculate_enu_to_ecef_rotation_matrix(
        qm.header.center,
        &quantized_mesh::Ellipsoid::wgs84(),
    )
    .transpose();
    let dist_enu =
        to_enu_matrix.transform_vector(&(qm.header.bounding_sphere.center - qm.header.center));
    println!(
        "Tile Centre -> Bounding Sphere Centre (ENU from tile centre): {:#?}",
        dist_enu
    );

    // Print the first triangle
    println!("0: {:#?}", qm.get_triangle(0, &bbox));

    Ok(())
}
