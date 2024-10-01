
use binrw::{prelude::*, BinRead, BinResult};
use binrw::helpers::until_eof;
use binrw::binread;
use geometry::{lerp, ECEFPoint3, Ellipsoid, GeodeticPoint2, GeodeticPoint3, GeodeticRectangle};

pub mod kml_writer;
pub mod svg_writer;
pub mod tile;
pub mod tiff_writer;
pub mod geometry;


// Constants
pub static UV_MAX_F64:f64 = 32767.0;
pub static UV_MAX_U16:u16 = 32767;


#[derive(Default)]
pub enum Projection {
    #[default]
    Epsg4326 = 1,
    Epsg3857 = 2,
}

#[derive(Debug, Default)]
pub enum TilingScheme {
    #[default]
    Tms = 1,
    Slippy = 2,
}



pub type TriangleVertexIndex = [usize; 3];



#[binread]
#[derive(Debug)]
#[br(little)]
pub struct BoundingSphere {
    #[br(parse_with = parse_point3)]
    pub center: ECEFPoint3<f64>,
    pub radius: f64,
}

/// Convert WorldCRS84Quad tile to bounding box (lat/lon).
/// Returns a LL/UR bounding box
pub fn tile_to_bbox_crs84(x: u32, y: u32, z: u32) -> GeodeticRectangle<f64>{
    let tiles_per_side = 2 << z; // 2 tiles at 0 level for WGS84
    let tile_size_deg = 360.0 / tiles_per_side as f64; // Tile size in degrees

    // Longitude bounds
    let min_lon = x as f64 * tile_size_deg - 180.0;
    let max_lon = (x as f64 + 1.0) * tile_size_deg - 180.0;

    // Latitude boundssc
    let min_lat = (y as f64) * tile_size_deg - 90.0;
    let max_lat = (y as f64 + 1.0) * tile_size_deg - 90.0;


    GeodeticRectangle {
        lower_left: GeodeticPoint2::from_degrees(min_lon, min_lat),
        upper_right: GeodeticPoint2::from_degrees(max_lon, max_lat),
    }
}

// Quantized Mesh Binary Handling

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct QuantizedMeshHeader {
    #[br(parse_with = parse_point3)]
    pub center: ECEFPoint3<f64>,

    pub min_height: f32,

    pub max_height: f32,

    pub bounding_sphere: BoundingSphere,

    #[br(parse_with = parse_point3)]
    pub horizon_occlusion_point: ECEFPoint3<f64>,
}

#[binread]
#[derive(Debug)]
#[br(little)]
    pub struct VertexData {
        pub vertex_count: u32,        
        
        #[br(parse_with = vertex_vec_decode, args (vertex_count))]
        pub u: Vec<u16>,
        
        #[br(parse_with = vertex_vec_decode, args (vertex_count))]
        pub v: Vec<u16>,
    
        #[br(parse_with = vertex_vec_decode, args (vertex_count))]
        pub height: Vec<u16>,

        #[br(align_before = if vertex_count > 65536 { 4 } else { 2 })]

        pub triangle_count: u32,

        #[br(parse_with = triangle_vec_decode, args ( triangle_count, vertex_count > 65536))]
        pub triangle_index: Vec<TriangleVertexIndex>,

        #[br(args {
        long: vertex_count > 65536
        })]
        pub edge_indices: EdgeIndices,
    }


#[binread]
#[derive(Debug)]
#[br(little)]
pub struct QuantizedMesh {
    #[br(ignore)]
    pub bounding_rectangle: GeodeticRectangle<f64>,

    #[br(ignore)]
    pub ellipsoid: Ellipsoid,

    #[br(ignore)]
    pub tiling_scheme: TilingScheme,

    pub header: QuantizedMeshHeader,
    pub vertex_data: VertexData,

    #[br(parse_with = until_eof)]
    pub extensions: Vec<Extension>,

}

impl QuantizedMesh {
    pub fn new(
        header: QuantizedMeshHeader,
        vertex_data: VertexData,
        extensions: Vec<Extension>,
        bounding_rectangle: GeodeticRectangle<f64>,
        ellipsoid: Ellipsoid,
        tiling_scheme: TilingScheme,
    ) -> Self {
        QuantizedMesh {
            header,
            vertex_data,
            extensions,
            bounding_rectangle,
            ellipsoid,
            tiling_scheme,
        }
    }

    pub fn vertex_as_geodetic_point3(
        &self,
        vertex_index: usize,
    ) -> GeodeticPoint3<f64> {
        let u_value = self.vertex_data.u[vertex_index] as f64;
        let v_value = self.vertex_data.v[vertex_index] as f64;
        let height_value = self.vertex_data.height[vertex_index] as f64;

        let ll = &self.bounding_rectangle.lower_left;
        let ur  = &self.bounding_rectangle.upper_right;

        let lat = lerp(ll.lon(), ur.lon(), &(u_value / UV_MAX_F64));
        let lon = lerp(ll.lat(), ur.lat(), &(v_value / UV_MAX_F64));
        let alt: f64 = lerp(&(self.header.min_height as f64), &(self.header.max_height as f64), &(height_value / UV_MAX_F64));

        GeodeticPoint3::new(lat, lon, alt)
    }

    pub fn vertices_as_geodetic_point3(
        &self,
    ) -> Vec<GeodeticPoint3<f64>> {
        let mut geodetic_coords = Vec::with_capacity(self.vertex_data.vertex_count as usize);
        for i in 0..self.vertex_data.vertex_count as usize {
            let lat_lon_height =
                self.vertex_as_geodetic_point3(i);
            geodetic_coords.push(lat_lon_height);
        }
        geodetic_coords
    }

    pub fn interpolated_height_vertices(
        &self,
    ) -> Vec<f64> {
        let mut heights: Vec<f64> = Vec::with_capacity(self.vertex_data.vertex_count as usize);
        for i in 0..self.vertex_data.vertex_count as usize {
            heights.push(lerp(&(self.header.min_height as f64), &(self.header.max_height as f64), &(self.vertex_data.height[i] as f64 / UV_MAX_F64)));
        }
        heights
    }
    
}



#[binread   ]
#[derive(Debug)]
#[br(little)]
#[derive(Default)]
#[br(import{ 
    long: bool}
)]
pub struct EdgeIndices{
    pub west_vertex_count: u32,
    #[br(parse_with = edge_index_decode, args(west_vertex_count, long))]
    pub west_indices: Vec<usize>,

    pub south_vertex_count: u32,
    #[br(parse_with = edge_index_decode, args(south_vertex_count, long))]
    pub south_indices: Vec<usize>,

    pub east_vertex_count: u32,
    #[br(parse_with = edge_index_decode, args(east_vertex_count, long))]
    pub east_indices: Vec<usize>,

    pub north_vertex_count: u32,
    #[br(parse_with = edge_index_decode, args(north_vertex_count, long))]
    pub north_indices: Vec<usize>,
}

#[binread]
#[derive(Debug)]
#[repr(u8)]
#[br(repr = u8)]
pub enum ExtensionID {
    TerrainLighting = 1,
    WaterMask = 2,
    Metadata = 4,
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct Extension {
    pub extension_id: u8,

    pub extension_length: u32, // bytes

    #[br(count = extension_length)]
    pub extension_data: Vec<u8>,
}

#[binrw::parser(reader, endian)]
fn parse_point3() -> BinResult<ECEFPoint3<f64>> {
    Ok(ECEFPoint3::new(
        <_>::read_options(reader, endian, ())?,
        <_>::read_options(reader, endian, ())?,
        <_>::read_options(reader, endian, ())?,
    ))
}

// An edge index can be vector of u16 or u32. Read either and return as vec<u32>
#[binrw::parser(reader,endian)]
pub fn edge_index_decode(count: u32, long: bool) -> binrw::BinResult<Vec<usize>> {
    let mut result = Vec::with_capacity(count as usize);
    
    // Read 'count' number of u16 values and convert to u32
    for _ in 0..count {
        if long {
            result.push(u32::read_options(reader, endian, ())? as usize);
        }
        else {
            result.push( u16::read_options(reader, endian, ())? as usize);
        }
    }
    
    Ok(result)
}


// Vertices are delta and zig-zag encoded as a vec<u16>
// This encoding can produce a negative index value
// In this case we know min u,v=0
#[binrw::parser(reader,endian)]
pub fn vertex_vec_decode(count: u32) ->  binrw::BinResult<Vec<u16>> {
    let mut res: Vec<u16> = Vec::with_capacity(count as usize);
    let mut val: u16 = 0;

    for _ in 0..count {
        
        let n = u16::read_options(reader, endian, ())?;
        let decoded = (n >>1 )as i16 ^ -((n  & 1) as i16);
        // Attempt to add decoded to val, checking for overflow
        
    // Attempt to add decoded to val, checking for overflow
    val = val
    .checked_add_signed(decoded)
    .ok_or_else(|| {
        binrw::Error::AssertFail {
            pos: reader.stream_position().unwrap(), // Include the stream position in the error
            message: format!("Overflow occurred while decoding value: {}", decoded),
        }
    })?;
        
       
        res.push(val);
    }
   Ok(res)
} 

// Triangle index is high watermark encoded
// Sanity check included as this encoding can 
#[binrw::parser(reader, endian)]
pub fn triangle_vec_decode(count: u32, long: bool) -> binrw::BinResult<Vec<[usize; 3]>> {
    let mut res: Vec<[usize; 3]> = Vec::with_capacity(count as usize);
    let mut highest: usize = 0;

    for _ in 0..count {
        // Initialize an array for the triangle
        let mut tri: TriangleVertexIndex = [0; 3];

        for j in 0..3 {
            let n: usize = if long {
                u32::read_options(reader, endian, ())? as usize
            } else {
                u16::read_options(reader, endian, ())? as usize
            };

            // bounds check
            tri[j] = highest.checked_sub(n)
            .ok_or_else(|| {
                binrw::Error::AssertFail {
                    pos: reader.stream_position().unwrap(), // Include the stream position in the error
                    message: format!(
                        "Invalid high watermark index encoding - highest: {} < value: {}",
                        highest, n
                    ),
                }
            })?;

            if n == 0 {
                highest += 1;
            }
        }

        res.push(tri);
    }

    Ok(res)
}
