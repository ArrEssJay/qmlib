use binrw::BinRead;
use std::path::{Path, PathBuf};

use crate::{
    geometry::{
        lerp, Ellipsoid, GeodeticPoint2, GeodeticPoint3, GeodeticRectangle
    },
    QuantizedMesh, UV_SIZE_F64,
};

use std::fs::File;

// Emums
#[derive(Default, Debug, PartialEq)]
pub enum CRS {
    #[default]
    Epsg4326,
    Epsg3857,
}

#[derive(Debug, Default)]
pub enum TilingScheme {
    #[default]
    Tms,
    Slippy,
}

pub struct QuantizedMeshTile {
    pub quantized_mesh: QuantizedMesh,
    pub zoom: u32,
    pub x: u32,
    pub y: u32,
    pub bounding_rectangle: GeodeticRectangle<f64>,
    pub ellipsoid: Ellipsoid,
    pub tiling_scheme: TilingScheme,
    pub crs: CRS,
}

impl QuantizedMeshTile {
    pub fn new(
        quantized_mesh: QuantizedMesh,
        zoom: u32,
        x: u32,
        y: u32,
        ellipsoid: Ellipsoid,
        tiling_scheme: TilingScheme,
        crs: CRS,
    ) -> Self {
        Self {
            quantized_mesh,
            zoom,
            x,
            y,
            bounding_rectangle: tile_to_bounding_rectangle(&crs, x, y, zoom).unwrap(),
            ellipsoid,
            tiling_scheme,
            crs,
        }
    }

    pub fn vertex_as_geodetic_point3(&self, vertex_index: usize) -> GeodeticPoint3<f64> {
        let u_value = f64::from(self.quantized_mesh.vertex_data.u[vertex_index]);
        let v_value = f64::from(self.quantized_mesh.vertex_data.v[vertex_index]);
        let height_value = f64::from(self.quantized_mesh.vertex_data.height[vertex_index]);

        let ll = &self.bounding_rectangle.lower_left;
        let ur = &self.bounding_rectangle.upper_right;

        let lat = lerp(ll.lon(), ur.lon(), &(u_value / (UV_SIZE_F64 -1.0)));
        let lon = lerp(ll.lat(), ur.lat(), &(v_value /  (UV_SIZE_F64 -1.0)));
        let alt: f64 = lerp(
            &f64::from(self.quantized_mesh.header.min_height),
            &f64::from(self.quantized_mesh.header.max_height),
            &(height_value /  (UV_SIZE_F64 -1.0)),
        );

        GeodeticPoint3::new(lat, lon, alt)
    }

    pub fn vertices_as_geodetic_point3(&self) -> Vec<GeodeticPoint3<f64>> {
        let mut geodetic_coords =
            Vec::with_capacity(self.quantized_mesh.vertex_data.vertex_count as usize);
        for i in 0..self.quantized_mesh.vertex_data.vertex_count as usize {
            let lat_lon_height = self.vertex_as_geodetic_point3(i);
            geodetic_coords.push(lat_lon_height);
        }
        geodetic_coords
    }
}

/// Convert `WorldCRS84Quad` tile to bounding box (lat/lon).
/// Returns a LL/UR bounding box
pub fn tile_to_bounding_rectangle(
    crs: &CRS,
    x: u32,
    y: u32,
    zoom: u32,
) -> Result<GeodeticRectangle<f64>, String> {
    if crs == &CRS::Epsg4326 {
        let tiles_per_side = 2 << zoom; // 2 tiles at 0 level for WGS84
        let tile_size_deg = 360.0 / tiles_per_side as f64; // Tile size in degrees

        // Longitude bounds
        let min_lon = x as f64 * tile_size_deg - 180.0;
        let max_lon = (x as f64 + 1.0) * tile_size_deg - 180.0;

        // Latitude boundssc
        let min_lat = (y as f64) * tile_size_deg - 90.0;
        let max_lat = (y as f64 + 1.0) * tile_size_deg - 90.0;

        Ok(GeodeticRectangle {
            lower_left: GeodeticPoint2::from_degrees(min_lon, min_lat),
            upper_right: GeodeticPoint2::from_degrees(max_lon, max_lat),
        })
    } else {
        Err(format!("Unsupported CRS {:?}", crs))
    }
}

pub fn decode_tile_from_path(path: &Path) -> Result<(usize, usize, usize), String> {
    // Extract the components from the path
    let components: Vec<&str> = path.iter().filter_map(|p| p.to_str()).collect();

    if components.len() < 3 {
        return Err("Path must contain at least three components for zoom, x, and y.".to_string());
    }

    // The last component should be the filename, which we will split to get the y value
    let file_name = components.last().ok_or("No file name found.")?;
    let y_str = file_name
        .split('.')
        .next()
        .ok_or("No y value found in the file name.")?;

    // Extract zoom, x from the path components
    let zoom = components[components.len() - 3]
        .parse::<usize>()
        .map_err(|_| "Failed to parse zoom value.")?;
    let x = components[components.len() - 2]
        .parse::<usize>()
        .map_err(|_| "Failed to parse x value.")?;
    let y = y_str
        .parse::<usize>()
        .map_err(|_| "Failed to parse y value.")?;

    Ok((zoom, x, y))
}

// TODO - Parse layer JSON and handle tile CRS parameters
pub fn load_quantized_mesh_tile(path: &PathBuf) -> Result<QuantizedMeshTile, String> {
    // Decode zoom, x, y from the file path
    let (decoded_zoom, decoded_x, decoded_y) = decode_tile_from_path(path)?;

    let zoom = decoded_zoom as u32;
    let x = decoded_x as u32;
    let y = decoded_y as u32;

    println!("Decoded tile params: zoom = {zoom}, x = {x}, y = {y}");

    // Open the file
    let mut file =
        File::open(path).map_err(|e| format!("Failed to open file {}: {:?}", path.display(), e))?;

    // Read the QuantizedMesh
    let qm: QuantizedMesh = QuantizedMesh::read_le(&mut file)
        .map_err(|e| format!("Failed to read quantized mesh: {e:?}"))?;

    // TODO: Determine georeferencing parameters from layer.json
    let qmtile = QuantizedMeshTile::new(
        qm,
        zoom,
        x,
        y,
        Ellipsoid::wgs84(),
        TilingScheme::Slippy,
        CRS::Epsg4326,
    );

    Ok(qmtile)
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;
    use crate::quantized_mesh_tile;

    #[test]
    fn test_load_tile() {
        // Assume you have a valid test tile path
        let path = PathBuf::from(format!("{}/test/terrain_data/a/15/59489/9692.terrain", env!("CARGO_MANIFEST_DIR")));       
        let result = quantized_mesh_tile::load_quantized_mesh_tile(&path);
        assert!(result.is_ok(), "Failed to load tile: {:?}", result.err());
    }
}