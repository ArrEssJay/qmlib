# Rust quantized-mesh-1.0 terrain format tools

This is a very rough and ready implementation of the [Cesium quantized-mesh-1.0](https://github.com/CesiumGS/quantized-mesh) terrain format in Rust, primarily written as a Rust learning exercise for myself.

At present binary parsing is (partially) implemented using `binrw`. A single tile can be unpacked and the uvh co-ordinates mapped from the ENU local tangent plane to a geodetic/projected frame. Only CRS84 (unprojected lat/lon) is currently handled.

## Example Applications
- `cargo run --bin tile_info <z/y/x.terrain>` print out tile header and georeferencing data. Expects tile to be located in a folder reflecting TMS hierachy.
- `cargo run --bin tile_export <z/y/x.terrain>` outputs SVG and KML files with vertices & mesh triangles.

## ToDo
- Write to a well supported mesh format/link to MDAL
- Read layer json and merge multiple tiles at edges
- Writing
- EPSG:3857 (Web Mercator) support
