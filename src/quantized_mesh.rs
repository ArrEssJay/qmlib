use binrw::helpers::until_eof;
use binrw::{binread, BinRead, BinResult};
use nalgebra::{Matrix4, Point3, Vector3};

/// Ellipsoid constants for WGS-84
#[derive(Debug)]
pub struct Ellipsoid {
    pub semi_major_axis: f64,
    pub flattening: f64,
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

    /// Calculate the surface normal at the given ECEF position (normalized vector)
    pub fn geodetic_surface_normal(&self, position: Point3<f64>) -> Vector3<f64> {
        let a = self.semi_major_axis;
        let e2 = self.eccentricity_squared();
        let (x, y, z) = (position.x, position.y, position.z);

        let p = (x * x + y * y).sqrt();
        let normal = Vector3::new(x / a, y / a, z / (a * (1.0 - e2).sqrt()));
        normal.normalize()
    }
}

/// Define a struct to hold the rectangle's vertices
#[derive(Debug)]
pub struct BoundingBox {
    lower_left: (f64, f64),  // (longitude, latitude)
    upper_right: (f64, f64), // (longitude, latitude)
}

impl BoundingBox {
    pub fn dimensions(&self) -> (f64, f64) {
        let width = self.upper_right.0 - self.lower_left.0; // Longitude difference
        let height = self.upper_right.1 - self.lower_left.1; // Latitude difference
        (width, height)
    }
}

/// Convert Cartesian ECEF coordinates to geodetic coordinates
pub fn ecef_to_geodetic(p: &Point3<f64>, el: &Ellipsoid) -> GeodeticPoint3 {
    let major = el.semi_major_axis;
    let minor = el.semi_minor_axis();

    let r = (p.x * p.x + p.y * p.y + p.z * p.z).sqrt();
    let e = (major * major - minor * minor).sqrt();
    let var = r * r - e * e;
    let u = (0.5 * var + 0.5 * (var * var + 4.0 * e * e * p.z * p.z).sqrt()).sqrt();

    let q = (p.x * p.x + p.y * p.y).sqrt();
    let mut beta = (u * u + e * e).sqrt() / u * p.z / q.atan();

    let eps = ((minor * u - major * beta.sqrt() + e * e) * beta.sin())
        / (major * beta.sqrt() / beta.cos() - e * e * beta.cos());
    beta += eps;

    let lat = (major / minor * beta.tan()).atan();
    let lon = p.y.atan2(p.x);
    let v1 = p.z - minor * beta.sin();
    let v2 = q - major * beta.cos();

    let alt =
        if (p.x * p.x / major / major) + (p.y * p.y / major / major) + (p.z * p.z / minor / minor)
            < 1.0
        {
            -(v1 * v1 + v2 * v2).sqrt()
        } else {
            (v1 * v1 + v2 * v2).sqrt()
        };

    GeodeticPoint3 {
        lat: lat.to_degrees(),
        lon: lon.to_degrees(),
        alt,
    }
}

/// Calculate the ENU to ECEF transformation matrix for a given ECEF reference point
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

#[binread]
#[derive(Debug)]
#[repr(u8)]
#[br(repr = u8)]
pub enum ExtensionID {
    TerrainLighting = 1,
    WaterMask = 2,
    Metadata = 4,
}

#[derive(Debug)]
pub struct GeodeticPoint3 {
    lat: f64, // Assuming LatLon is a f64
    lon: f64, // Assuming LatLon is a f64
    alt: f64, // Assuming Meter is a f64
}

pub trait ToDegrees {
    fn to_degrees(self) -> GeodeticPoint3;
}

impl ToDegrees for GeodeticPoint3 {
    fn to_degrees(self) -> GeodeticPoint3 {
        GeodeticPoint3 {
            lat: self.lat.to_degrees(),
            lon: self.lon.to_degrees(),
            alt: self.alt, // Altitude may not need conversion
        }
    }
}

pub trait ToRadians {
    fn to_radians(self) -> f64; // Assuming Radian is represented as f64
}

impl ToRadians for f64 {
    fn to_radians(self) -> f64 {
        self * std::f64::consts::PI / 180.0
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

pub trait EdgeIndices {}

impl EdgeIndices for EdgeIndicesShort {}

impl EdgeIndices for EdgeIndicesLong {}

#[binread]
#[derive(Debug)]
#[br(little)]
#[derive(Default)]
pub struct EdgeIndicesShort {
    #[br(dbg)]
    pub westVertexCount: u32,
    #[br(count = westVertexCount)]
    #[br(map = |vals: Vec<u16>| vals.into_iter().map(|v| v as u32).collect())]
    // Map u16 to u32
    pub westIndices: Vec<u32>,

    #[br(dbg)]
    pub southVertexCount: u32,
    #[br(count = southVertexCount)]
    #[br(map = |vals: Vec<u16>| vals.into_iter().map(|v| v as u32).collect())]
    pub southIndices: Vec<u32>,

    #[br(dbg)]
    pub eastVertexCount: u32,
    #[br(count = eastVertexCount)]
    #[br(map = |vals: Vec<u16>| vals.into_iter().map(|v| v as u32).collect())]
    pub eastIndices: Vec<u32>,

    #[br(dbg)]
    pub northVertexCount: u32,
    #[br(count = northVertexCount)]
    #[br(map = |vals: Vec<u16>| vals.into_iter().map(|v| v as u32).collect())]
    pub northIndices: Vec<u32>,
}

#[binread]
#[derive(Debug)]
#[br(little)]
#[derive(Default)]
pub struct EdgeIndicesLong {
    #[br(dbg)]
    pub westVertexCount: u32,
    #[br(count = westVertexCount)]
    pub westIndices: Vec<u32>,

    #[br(dbg)]
    pub southVertexCount: u32,
    #[br(count = southVertexCount)]
    pub southIndices: Vec<u32>,

    #[br(dbg)]
    pub eastVertexCount: u32,
    #[br(count = eastVertexCount)]
    pub eastIndices: Vec<u32>,

    #[br(dbg)]
    pub northVertexCount: u32,
    #[br(count = northVertexCount)]
    pub northIndices: Vec<u32>,
}

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct VertexData {
    #[br(dbg)]
    pub vertex_count: u32,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i) )]
    pub u: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i))]
    pub v: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i))]
    pub height: Vec<i32>,

    #[br(dbg)]
    pub triangleCount: u32,

    #[br(align_before = 4)]
    #[br(count = triangleCount * 3)]
    #[br(if(vertex_count > 65536) )]
    #[br(map = |i:  Vec<u32>| Some(i.high_watermark_decode()))]
    pub index_data_long: Option<Vec<usize>>,

    #[br(align_before = 2)]
    #[br(count = triangleCount * 3)]
    #[br(if(vertex_count <= 65536))]
    #[br(map = |i: Vec<u16>| {
        // Convert Vec<u16> to Vec<u32>
        let u32_vec: Vec<u32> = i.into_iter().map(|x| x as u32).collect();
        // Apply high_watermark_decode to Vec<u32> and wrap in Some
        Some(u32_vec.high_watermark_decode())
    })]
    pub index_data_short: Option<Vec<usize>>,

    #[br(align_before = 4)]
    #[br(if(vertex_count > 65536) )]
    pub edge_indices_long: Option<EdgeIndicesLong>,

    #[br(align_before = 2)]
    #[br(if(vertex_count <= 65536) )]
    pub edge_indices_short: Option<EdgeIndicesShort>,
}

trait HighWatermarkDecode {
    fn high_watermark_decode(&self) -> Vec<usize>;
}

impl HighWatermarkDecode for Vec<u32> {
    fn high_watermark_decode(&self) -> Vec<usize> {
        let mut decoded: Vec<usize> = Vec::with_capacity(self.len());
        let mut highest: usize = 0;

        for &i in self.iter() {
            let i_as_usize = i as usize;
            decoded.push(highest - i_as_usize);
            if i == 0 {
                highest += 1;
            }
        }
        decoded
    }
}

fn vec_zigzag_decode(src: Vec<u16>) -> Vec<i32> {
    let mut res: Vec<i32> = Vec::with_capacity(src.len());
    let mut n: i32 = 0;
    for i in src.iter() {
        let j: i32 = *i as i32;
        // zigzag decode
        n += (j >> 1) ^ (-(j & 1));
        res.push(n);
    }
    res
}

#[binread]
#[derive(Debug)]
#[br(little)]
pub struct QuantizedMeshHeader {
    #[br(dbg)]
    #[br(parse_with = parse_point3)]
    pub center: Point3<f64>,

    #[br(dbg)]
    pub min_height: f32,

    #[br(dbg)]
    pub max_height: f32,

    #[br(dbg)]
    pub bounding_sphere: BoundingSphere,

    #[br(dbg)]
    #[br(parse_with = parse_point3)]
    pub horizon_occlusion_point: Point3<f64>,
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct Extension {
    #[br(dbg)]
    pub extensionId: u8,
    #[br(dbg)]
    pub extensionLength: u32, // bytes

    #[br(count = extensionLength)]
    pub extensionData: Vec<u8>,
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct QuantizedMesh {
    pub header: QuantizedMeshHeader,
    pub vertex_data: VertexData,

    #[br(parse_with = until_eof)]
    pub extensions: Vec<Extension>,
}

impl QuantizedMesh {
    fn get_interpolated_vertex(&self, vertex_index: usize, bbox: &BoundingBox) -> Point3<f64> {
        let u_value = self.vertex_data.u[vertex_index] as f64;
        let v_value = self.vertex_data.v[vertex_index] as f64;
        let height_value = self.vertex_data.height[vertex_index] as f64;

        let (min_lon, min_lat) = bbox.lower_left; // Lower-left corner
        let (max_lon, max_lat) = bbox.upper_right; // Upper-right corner

        // Lerp for longitude (x) and latitude (y)
        let x = lerp(min_lon, max_lon, u_value / 32767.0);
        let y = lerp(min_lat, max_lat, v_value / 32767.0);

        // Lerp for height (z)
        let z = lerp(
            self.header.min_height as f64,
            self.header.max_height as f64,
            height_value / 32767.0,
        );

        Point3::new(x, y, z)
    }
    pub fn get_triangle(
        &self,
        triangle_index: usize,
        bbox: &BoundingBox,
    ) -> Option<[Point3<f64>; 3]> {
        // long or short index?
        let index_data = match (
            self.vertex_data.index_data_long.as_ref(),
            self.vertex_data.index_data_short.as_ref(),
        ) {
            (Some(long), _) => long,
            (None, Some(short)) => short,
            (None, None) => return None, // Both options are None
        };

        // check triangle index bounds
        //assert!(triangle_index * 3 + 2 >= index_data.len());

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

        // Calculate interpolated coordinates for each vertex
        let triangle: [nalgebra::OPoint<f64, nalgebra::Const<3>>; 3] = [
            self.get_interpolated_vertex(triangle_vertices[0], &bbox),
            self.get_interpolated_vertex(triangle_vertices[1], &bbox),
            self.get_interpolated_vertex(triangle_vertices[2], &bbox),
        ];

        Some(triangle)
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

// Linear interpolation function
fn lerp(min_value: f64, max_value: f64, t: f64) -> f64 {
    min_value + t * (max_value - min_value)
}

/// Convert WorldCRS84Quad tile to bounding box (lat/lon).
/// Returns the bounding box (min_lon, min_lat, max_lon, max_lat).
pub fn tile_to_bbox_crs84(x: u32, y: u32, z: u32) -> (BoundingBox) {
    let tiles_per_side = 2 << z; // Number of tiles along one side at this zoom level
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
