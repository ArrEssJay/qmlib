# Rust Quantized-Mesh 1.0 Terrain Format Tools

This is a native Rust implementation of tooling for the [Cesium Quantized-Mesh 1.0](https://github.com/CesiumGS/quantized-mesh) terrain format. It is and shall likely remain incomplete, primarily undertaken by myself as a Rust learning exercise.

## Features
- **Binary Reading**: Partial reading of a single mesh tile (extensions are not yet handled, nor is writing).
- **Georeferencing**: Converts local UVH coordinates to Geodetic LLH coordinates (WGS84 only).
- **Exporting**: Outputs vector paths to SVG and KML formats.
- **Multithreaded Rasterisation**: Interpolates elevation over triangle faces to create a DEM (GeoTIFF).

Currently, reading of `layer.json` and the handling of gridded tile sets are not implemented. The tile loader enforces a WGS84 Geographic CRS with a TMS tiling structure for any loaded tile.

Georeferencing parameters are derived from the `{z}`, `{y}`, `{x}.terrain` components of the supplied tile path in the tile loader.

For performance, georeferencing and rasterisation use `nalgebra` matrix methods, and the raster pipeline is multithreaded.

## Example Applications
- `cargo run --bin tile_info <z/y/x.terrain>`: Prints out tile header and georeferencing data. Expects the tile to be located in a folder reflecting the TMS hierarchy.
- `cargo run --bin export_svg <z/y/x.terrain>`: Outputs SVG with unreferenced UV vertices and triangle edges.
- `cargo run --bin export_kml <z/y/x.terrain>`: Outputs KML with geodetic latitude, longitude, and height vertices and triangle edges.
- `cargo run --bin export_geotiff <z/y/x.terrain> <scale shift>`: Outputs a GeoTIFF. Output resolution is determined by the scale shift parameter, right shifting vertex UV co-ordinates N bits. Spatial resolution will depend on the area covered by the tile.

## Tests/Benchmarks
- `cargo test` - Run tests
- `cargo bench --release` - Run Criterion benchmarks with release optimisation flags

## To Do
- Write to a well-supported mesh format or link to MDAL.
- Read `layer.json` and handling of tiled layers in general.
- Support extensions and writing.
- Implement EPSG:3857 (Web Mercator) support.
- Test coverage.
- Performance improvements.
