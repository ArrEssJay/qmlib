# Rust quantized-mesh-1.0 terrain format tools

This is a very rough and ready implementation of the Cesium  quantized-mesh-1.0 terrain format in Rust, primarily written as a Rust learning exercise for myself.

At present only parsing is (partially) implemented. A single tile can be unpacked and the uvh co-ordinates mapped from the ENU local tangent plane to a geodetic/projected frame. Only CRS84 (unprojected lat/lon) is currently handled.

The `dump_tri` application demonstrates usage of the library.