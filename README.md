# Rust Quantized-Mesh 1.0 Terrain Format Tools

This is a native Rust implementation of tooling for the [Cesium Quantized-Mesh 1.0](https://github.com/CesiumGS/quantized-mesh) terrain format. It is and shall likely remain incomplete, primarily undertaken by myself as a Rust learning exercise.

## Features
- **Binary Reading**: Partial reading of a single mesh tile (extensions are not yet handled, nor is writing).
- **Georeferencing**: Converts local UVH coordinates to Geodetic LLH coordinates (WGS84 only).
- **Exporting**: Outputs vector paths to SVG and KML formats.
- **Rasterization**: Interpolates triangle faces to create a DEM (GeoTIFF).

Currently, reading of `layer.json` and the handling of gridded tile sets are not implemented. The tile loader enforces a WGS84 Geographic CRS with a TMS tiling structure for any loaded tile.

Georeferencing parameters are derived from the `{z}`, `{y}`, `{x}.terrain` components of the supplied tile path in the tile loader.

For performance, georeferencing and rasterization utilize `nalgebra` matrix methods, though there is likely room for optimization.

## Example Applications
- `cargo run --bin tile_info <z/y/x.terrain>`: Prints out tile header and georeferencing data. Expects the tile to be located in a folder reflecting the TMS hierarchy.
- `cargo run --bin export_svg <z/y/x.terrain>`: Outputs SVG with unreferenced UV vertices and triangle edges.
- `cargo run --bin export_kml <z/y/x.terrain>`: Outputs KML with geodetic latitude, longitude, and height vertices and triangle edges.
- `cargo run --bin export_geotiff <z/y/x.terrain> <scale shift>`: Outputs a GeoTIFF. Output resolution is determined by the scale shift parameter, right shifting vertex UV co-ordinates N bits. Spatial resolution will depend on the area covered by the tile.

For the WGS84 Geographic case, approximate (equatorial) raster cell sizes (metres) at common TMS levels for a given output scale are shown here:

| N |   X,Y Resolution | **TMS Levels** |              |              |              |              |              |              |              |              |              |               |               |               |               |               |               |               |               |               |               |
|---|------------------|----------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|--------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|
|   |                  | Level 0        | Level 1      | Level 2      | Level 3      | Level 4      | Level 5      | Level 6      | Level 7      | Level 8      | Level 9      | Level 10      | Level 11      | Level 12      | Level 13      | Level 14      | Level 15      | Level 16      | Level 17      | Level 18      | Level 19      | Level 20      |
| 0 | 32768            | 156412.0       | 78206.0      | 39103.0      | 19551.5      | 9775.8       | 4887.9       | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7         | 76.4          | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 1 | 16384            | 78206.0        | 39103.0      | 19551.5      | 9775.8       | 4887.9       | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4          | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 2 | 8192             | 39103.0        | 19551.5      | 9775.8       | 4887.9       | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4         | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 3 | 4096             | 19551.5        | 9775.8       | 4887.9       | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4         | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 4 | 2048             | 9775.8         | 4887.9       | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4         | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 5 | 1024             | 4887.9         | 2443.9       | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4         | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |
| 6 | 512              | 2443.9         | 1221.9       | 610.9        | 305.4        | 152.7        | 76.4         | 38.2          | 19.1          | 9.5           | 4.8           | 2.4           | 1.2           | 0.6           | 0.3           | 0.15          |

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
