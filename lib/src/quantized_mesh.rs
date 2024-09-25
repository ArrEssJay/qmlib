use binrw::{prelude::*, BinRead, BinResult};
use binrw::helpers::until_eof;
use binrw::binread;
use nalgebra::{Matrix4, Point3, Vector3};

pub mod kml_writer;
pub mod svg_writer;
pub mod tile;
/// Ellipsoids
#[derive(Debug, Default)]
pub struct Ellipsoid {
    pub semi_major_axis: f64, // Semi-major axis (equatorial radius) in meters
    pub flattening: f64,      // Semi-minor axis (polar radius) in meters
}

impl Ellipsoid {
    /// Create a WGS-84 Ellipsoid
    pub fn wgs84() -> Self {
        Self {
            semi_major_axis: 6378137.0, // meters
            flattening: 1.0 / 298.257223563,
        }
    }

    /// Compute the first eccentricity squared (e²)
    pub fn eccentricity_squared(&self) -> f64 {
        self.flattening * (2.0 - self.flattening)
    }

    /// Compute the semi-minor axis (b)
    pub fn semi_minor_axis(&self) -> f64 {
        self.semi_major_axis * (1.0 - self.flattening)
    }

    /// Calculate the surface normal vector at the given ECEF position
    pub fn geodetic_surface_normal(&self, p: Point3<f64>) -> Vector3<f64> {
        let normal = Vector3::new(
            p.x / (self.semi_major_axis * self.semi_major_axis),
            p.y / (self.semi_major_axis * self.semi_major_axis),
            p.z / (self.semi_minor_axis() * (self.semi_minor_axis())),
        );
        // Normalize the vector to get the unit normal
        normal.normalize()
    }
}

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

#[derive(Debug, Default)]
pub struct BoundingBox {
    lower_left: (f64, f64),
    upper_right: (f64, f64),
}

impl BoundingBox {
    pub fn dimensions(&self) -> (f64, f64) {
        let width = self.upper_right.0 - self.lower_left.0;
        let height = self.upper_right.1 - self.lower_left.1;
        (width, height)
    }
}

pub type CartesianPoint3 = Point3<f64>;
pub type GeodeticPoint3 = Point3<f64>;

pub trait ToDegrees {
    fn to_degrees(self) -> GeodeticPoint3;
}

impl ToDegrees for GeodeticPoint3 {
    fn to_degrees(self) -> GeodeticPoint3 {
        GeodeticPoint3::new(self[0].to_degrees(), self[1].to_degrees(), self[2])
    }
}

pub trait ToRadians {
    fn to_radians(self) -> GeodeticPoint3;
}

impl ToRadians for GeodeticPoint3 {
    fn to_radians(self) -> GeodeticPoint3 {
        GeodeticPoint3::new(self[0].to_radians(), self[1].to_radians(), self[2])
    }
}

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct BoundingSphere {
    #[br(parse_with = parse_point3)]
    pub center: Point3<f64>,
    pub radius: f64,
}

/// Convert Cartesian ECEF coordinates to geodetic coordinates
/// Rey-Jer You non-iterative method
pub fn ecef_to_geodetic(p: &Point3<f64>, el: &Ellipsoid) -> GeodeticPoint3 {
    let major = el.semi_major_axis;
    let minor = el.semi_minor_axis();

    let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    let e = (major * major - minor * minor).sqrt();
    let var = r * r - e * e;
    let u = (0.5 * var + 0.5 * (var * var + 4.0 * e * e * p.z * p.z).sqrt()).sqrt();

    let q = (p.x * p.x + p.y * p.y).sqrt();
    let hu_e = (u * u + e * e).sqrt();
    let mut beta = (hu_e / u * p.z / q).atan();

    let eps = ((minor * u - major * hu_e + e * e) * beta.sin())
        / (major * hu_e / beta.cos() - e * e * beta.cos());
    beta += eps;

    let lat = (major / minor * beta.tan()).atan();
    let lon = p.y.atan2(p.x);

    let v1 = p.z - minor * beta.sin();
    let v2 = q - major * beta.cos();
    let alt;

    let inside =
        (p.x * p.x / major / major) + (p.y * p.y / major / major) + (p.z * p.z / minor / minor)
            < 1.0;
    if inside {
        alt = -(v1 * v1 + v2 * v2).sqrt();
    } else {
        alt = (v1 * v1 + v2 * v2).sqrt();
    };

    GeodeticPoint3::new(lat, lon, alt)
}

/// Calculate the 4x4 ENU to ECEF rotation matrix for a given ECEF reference point
pub fn calculate_enu_to_ecef_rotation_matrix(
    ecef_position: Point3<f64>,
    ellipsoid: &Ellipsoid,
) -> Matrix4<f64> {
    let up = ellipsoid.geodetic_surface_normal(ecef_position);
    let mut east = Vector3::new(-ecef_position.y, ecef_position.x, 0.0);
    east.normalize_mut();
    let north = up.cross(&east).normalize();

    Matrix4::new(
        east.x,
        north.x,
        up.x,
        0.0,
        east.y,
        north.y,
        up.y,
        0.0,
        east.z,
        north.z,
        up.z,
        0.0,
        ecef_position.x,
        ecef_position.y,
        ecef_position.z,
        1.0,
    )
}

// Linear interpolation function
fn lerp(min_value: f64, max_value: f64, t: f64) -> f64 {
    min_value + t * (max_value - min_value)
}

/// Convert WorldCRS84Quad tile to bounding box (lat/lon).
/// Returns a LL/UR bounding box
pub fn tile_to_bbox_crs84(x: u32, y: u32, z: u32) -> BoundingBox {
    let tiles_per_side = 2 << z; // 2 tiles at 0 level for WGS84
    let tile_size_deg = 360.0 / tiles_per_side as f64; // Tile size in degrees

    // Longitude bounds
    let min_lon = x as f64 * tile_size_deg - 180.0;
    let max_lon = (x as f64 + 1.0) * tile_size_deg - 180.0;

    // Latitude boundss
    let min_lat = (y as f64) * tile_size_deg - 90.0;
    let max_lat = (y as f64 + 1.0) * tile_size_deg - 90.0;

    BoundingBox {
        lower_left: (min_lon, min_lat),
        upper_right: (max_lon, max_lat),
    }
}

// Quantized Mesh Binary Handling

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct QuantizedMeshHeader {
    #[br(parse_with = parse_point3)]
    pub center: CartesianPoint3,

    pub min_height: f32,

    pub max_height: f32,

    pub bounding_sphere: BoundingSphere,

    #[br(parse_with = parse_point3)]
    pub horizon_occlusion_point: CartesianPoint3,
}

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct VertexData {
    pub vertex_count: u32,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| i.vertex_vec_decode() )]
    pub u: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| i.vertex_vec_decode())]
    pub v: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| i.vertex_vec_decode())]
    pub height: Vec<i32>,

    pub triangle_count: u32,

    #[br(align_before = 4)]
    #[br(count = triangle_count * 3)]
    #[br(if(vertex_count > 65536) )]
    #[br(map = |i:  Vec<u32>| Some(i.high_watermark_decode()))]
    pub index_data_long: Option<Vec<usize>>,

    #[br(align_before = 2)]
    #[br(count = triangle_count * 3)]
    #[br(if(vertex_count <= 65536))]
    #[br(map = |i: Vec<u16>| {
        let u32_vec: Vec<u32> = i.into_iter().map(|x| x as u32).collect();
        Some(u32_vec.high_watermark_decode())
    })]
    pub index_data_short: Option<Vec<usize>>,

    #[br(align_before = if vertex_count > 65536 { 4 } else { 2 })]
        #[br(args {
        long: vertex_count > 65536
    })]
    pub edge_indices: EdgeIndices,
}

impl VertexData {
    pub fn index_data(&self) -> Option<&[usize]> {
        // Check if long or short index data exists
        let index_data = match (
            self.index_data_long.as_ref(),
            self.index_data_short.as_ref(),
        ) {
            (Some(long), _) => long.as_slice(),
            (None, Some(short)) => short.as_slice(),
            (None, None) => return None,
        };

        Some(index_data)
    }

    pub fn to_geodetic(
        &self,
        bounding_box: &BoundingBox,
        header: &QuantizedMeshHeader,
    ) -> Vec<[f64; 3]> {
        let mut geodetic_coords = Vec::with_capacity(self.vertex_count as usize);
        for i in 0..self.vertex_count as usize {
            let lat_lon_height =
                self.to_geodetic_vertex(i, bounding_box, header.min_height, header.max_height);
            geodetic_coords.push(lat_lon_height);
        }
        geodetic_coords
    }

    pub fn to_geodetic_vertex(
        &self,
        vertex_index: usize,
        bbox: &BoundingBox,
        min_height: f32,
        max_height: f32,
    ) -> [f64; 3] {
        let u_value = self.u[vertex_index] as f64;
        let v_value = self.v[vertex_index] as f64;
        let height_value = self.height[vertex_index] as f64;

        let (min_lon, min_lat) = bbox.lower_left;
        let (max_lon, max_lat) = bbox.upper_right;

        let x = lerp(min_lon, max_lon, u_value / 32767.0);
        let y = lerp(min_lat, max_lat, v_value / 32767.0);
        let z = lerp(min_height as f64, max_height as f64, height_value / 32767.0);

        [x, y, z]
    }
}

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct QuantizedMesh {
    #[br(ignore)]
    pub bounding_box: BoundingBox,

    #[br(ignore)]
    pub ellipsoid: Ellipsoid,

    #[br(ignore)]
    pub tiling_scheme: TilingScheme,

    pub header: QuantizedMeshHeader,
    pub vertex_data: VertexData,

    #[br(parse_with = until_eof)]
    pub extensions: Vec<Extension>,

    #[br(ignore)]
    pub geodetic_vertices: Vec<[f64; 3]>,
}

impl QuantizedMesh {
    pub fn new(
        header: QuantizedMeshHeader,
        vertex_data: VertexData,
        extensions: Vec<Extension>,
        bounding_box: BoundingBox,
        ellipsoid: Ellipsoid,
        tiling_scheme: TilingScheme,
    ) -> Self {
        let geodetic_vertices = vertex_data.to_geodetic(&bounding_box, &header);
        QuantizedMesh {
            header,
            vertex_data,
            extensions,
            bounding_box,
            ellipsoid,
            tiling_scheme,
            geodetic_vertices,
        }
    }

    pub fn to_geodetic_vertex(&self, vertex_index: usize, bbox: &BoundingBox) -> [f64; 3] {
        self.vertex_data.to_geodetic_vertex(
            vertex_index,
            bbox,
            self.header.min_height,
            self.header.max_height,
        )
    }

    pub fn to_point3(point: [f64; 3]) -> Point3<f64> {
        Point3::new(point[0], point[1], point[2])
    }

    pub fn get_triangle(
        &self,
        triangle_index: usize,
        use_geodetic: bool,
    ) -> Option<[Point3<f64>; 3]> {
        let index_data = self.vertex_data.index_data()?;

        let triangle_vertices = [
            index_data[triangle_index * 3],
            index_data[triangle_index * 3 + 1],
            index_data[triangle_index * 3 + 2],
        ];

        // Check vertex index bounds
        if triangle_vertices
            .iter()
            .any(|&x| x >= self.vertex_data.vertex_count as usize)
        {
            return None; // Out of bounds
        }

        // Use geodetic or raw UV based on the flag
        let triangle: [Point3<f64>; 3] = if use_geodetic {
            [
                Self::to_point3(self.geodetic_vertices[triangle_vertices[0]]),
                Self::to_point3(self.geodetic_vertices[triangle_vertices[1]]),
                Self::to_point3(self.geodetic_vertices[triangle_vertices[2]]),
            ]
        } else {
            [
                Point3::new(
                    self.vertex_data.u[triangle_vertices[0]] as f64,
                    self.vertex_data.v[triangle_vertices[0]] as f64,
                    self.vertex_data.height[triangle_vertices[0]] as f64,
                ),
                Point3::new(
                    self.vertex_data.u[triangle_vertices[1]] as f64,
                    self.vertex_data.v[triangle_vertices[1]] as f64,
                    self.vertex_data.height[triangle_vertices[1]] as f64,
                ),
                Point3::new(
                    self.vertex_data.u[triangle_vertices[2]] as f64,
                    self.vertex_data.v[triangle_vertices[2]] as f64,
                    self.vertex_data.height[triangle_vertices[2]] as f64,
                ),
            ]
        };

        Some(triangle)
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
    #[br(dbg)]
    pub west_vertex_count: u32,
    #[br(dbg)]
    #[br(parse_with = read_u16_u32_as_u32, args(west_vertex_count as usize, long))]
    pub west_indices: Vec<u32>,

    #[br(dbg)]
    pub south_vertex_count: u32,
    #[br(parse_with = read_u16_u32_as_u32, args(south_vertex_count as usize, long))]
    pub south_indices: Vec<u32>,

    #[br(dbg)]
    pub east_vertex_count: u32,
    #[br(parse_with = read_u16_u32_as_u32, args(east_vertex_count as usize, long))]
    pub east_indices: Vec<u32>,

    #[br(dbg)]
    pub north_vertex_count: u32,
    #[br(dbg)]
    #[br(parse_with = read_u16_u32_as_u32, args(east_vertex_count as usize, long))]
    pub north_indices: Vec<u32>,
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

trait HighWatermarkDecode {
    fn high_watermark_decode(&self) -> Vec<usize>;
}

impl HighWatermarkDecode for Vec<u32> {
    fn high_watermark_decode(&self) -> Vec<usize> {
        let mut decoded: Vec<usize> = Vec::with_capacity(self.len());
        let mut highest: usize = 0;

        for &i in self.iter() {
            decoded.push(highest - (i as usize));
            if i == 0 {
                highest += 1;
            }
        }
        decoded
    }
}

pub fn zigzag_decode(n: i32) -> i32 {
    (n >> 1) ^ (-(n & 1))
}

trait VertexVecDecode {
    fn vertex_vec_decode(&self) -> Vec<i32>;
}
impl VertexVecDecode for Vec<u16> {
    fn vertex_vec_decode(&self) -> Vec<i32> {
        let mut res: Vec<i32> = Vec::with_capacity(self.len());
        let mut val = 0;

        for &i in self.iter() {
            let n = i as i32;
            val += zigzag_decode(n);
            res.push(val);
        }
        res
    }
}

#[binrw::parser(reader, endian)]
fn parse_point3() -> BinResult<Point3<f64>> {
    Ok(Point3::new(
        <_>::read_options(reader, endian, ())?,
        <_>::read_options(reader, endian, ())?,
        <_>::read_options(reader, endian, ())?,
    ))
}

// #[binrw::parser(reader, endian)]
// pub fn read_u16_as_u32() -> binrw::BinResult<u32> {
//     let val = u16::read_options(reader, endian, ())?;
//     println!("Read u16 bytes: {:?}", val);
//     Ok(val as u32) // Converts u16 to u32
// }


#[binrw::parser(reader,endian)]
pub fn read_u16_u32_as_u32(count: usize, long: bool) -> binrw::BinResult<Vec<u32>> {
    let mut result = Vec::with_capacity(count);
    
    // Read 'count' number of u16 values and convert to u32
    for _ in 0..count {
        if long {
            let val = u32::read_options(reader, endian, ())?;
            result.push(val);
        }
        else {
            let val = u16::read_options(reader, endian, ())?;
            result.push(val as u32);
        }
    }
    
    Ok(result)
}