# Rust quantized-mesh-1.0 terrain format tools

This is a very rough and ready implementation of the [Cesium quantized-mesh-1.0](https://github.com/CesiumGS/quantized-mesh) terrain format in Rust, primarily written as a Rust learning exercise for myself.

At present binary parsing is (partially) implemented using `binrw`. A single tile can be unpacked and the uvh co-ordinates mapped from the ENU local tangent plane to a geodetic/projected frame. Only CRS84 (unprojected lat/lon) is currently handled.

## Running
- `cargo run dump_tri` will build and print out some data from a supplied tile

## ToDo
- Write to a well supported mesh format/link to MDAL
- Read layer json and merge multiple tiles at edges
- Writing
- EPSG:3857 (Web Mercator) support
- Visualisation