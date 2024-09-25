use std::path::PathBuf;
use binrw::BinRead;

use crate::{tile_to_bbox_crs84, Ellipsoid, QuantizedMesh, TilingScheme};

use std::fs::File;


pub struct Tile {
    pub quantized_mesh: QuantizedMesh,
    pub zoom: u32,
    pub x: u32,
    pub y: u32,
}


pub fn decode_tile_from_path(path: &PathBuf) -> Result<(usize, usize, usize), String> {
    // Extract the components from the path
    let components: Vec<&str> = path
        .iter()
        .filter_map(|p| p.to_str())
        .collect();

    if components.len() < 3 {
        return Err("Path must contain at least three components for zoom, x, and y.".to_string());
    }

    // The last component should be the filename, which we will split to get the y value
    let file_name = components.last().ok_or("No file name found.")?;
    let y_str = file_name.split('.').next().ok_or("No y value found in the file name.")?;

    // Extract zoom, x from the path components
    let zoom = components[components.len() - 3].parse::<usize>()
        .map_err(|_| "Failed to parse zoom value.")?;
    let x = components[components.len() - 2].parse::<usize>()
        .map_err(|_| "Failed to parse x value.")?;
    let y = y_str.parse::<usize>()
        .map_err(|_| "Failed to parse y value.")?;

    Ok((zoom, x, y))
}

pub fn load_quantized_mesh(path: &PathBuf) -> Result<Tile, String> {
    // Decode zoom, x, y from the file path
    let (decoded_zoom, decoded_x, decoded_y) = decode_tile_from_path(path).map_err(|e| e)?;

    let zoom = decoded_zoom as u32;
    let x = decoded_x as u32;
    let y = decoded_y as u32;

    println!("Decoded tile params: zoom = {}, x = {}, y = {}", zoom, x, y);

    // Open the file
    let mut file = File::open(path).map_err(|e| {
        format!("Failed to open file {}: {:?}", path.display(), e)
    })?;

    // Read the QuantizedMesh
    let qmdata: QuantizedMesh = QuantizedMesh::read_le(&mut file).map_err(|e| {
        format!("Failed to read quantized mesh: {:?}", e)
    })?;

    let qm = QuantizedMesh::new(
        qmdata.header,
        qmdata.vertex_data,
        qmdata.extensions,
        tile_to_bbox_crs84(x, y, zoom),
        Ellipsoid::wgs84(),
        TilingScheme::Slippy,
    );

    Ok(Tile {
        quantized_mesh: qm,
        zoom,
        x,
        y,
    })
}
