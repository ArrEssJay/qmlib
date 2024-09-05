use std::ops::{Add, Sub};

use binrw::{
    binread,
    BinRead,
};
use itertools::{izip, TupleWindows, Tuples};
use map_3d::Ellipsoid;
pub type Meter<T> = T;
pub type Radian<T> = T;
pub type Degree<T> = T;
pub type Vertex3<T> = (T, T, T);

type Triangle<'a, T> = (&'a T, &'a T, &'a T);
//pub type Triangle<T> = (T, T, T);

#[derive(Debug)]
#[binread]
pub struct CartesianPoint3<Meter>
where
    Meter: for<'a> BinRead<Args<'a> = ()>, // Use an explicit lifetime 'a
                                           //     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
                                           //     T needs bounds to derive BinRead/BinWrite
{
    x: Meter,
    y: Meter,
    z: Meter,
}

impl CartesianPoint3<f64> {
    fn ecef_to_geodetic(self) -> GeodeticPoint3<f64, f64> {
        localecef2geodetic(self, map_3d::Ellipsoid::GRS80)
    }
}

pub struct GeodeticPoint3<LatLon, Meter> {
    lat: LatLon,
    lon: LatLon,
    alt: Meter,
}

pub trait ToDegrees<T> {
    fn to_degrees(self) -> Degree<T>;
}

impl ToDegrees<f64> for Radian<f64> {
    fn to_degrees(self) -> Degree<f64> {
        self.to_degrees()
    }
}

impl<LatLon, Meter> GeodeticPoint3<LatLon, Meter>
where
    LatLon: ToDegrees<f64>,
{
    fn to_degrees(self) -> GeodeticPoint3<Degree<f64>, Meter> {
        GeodeticPoint3 {
            lat: self.lat.to_degrees(),
            lon: self.lon.to_degrees(),
            alt: self.alt,
        }
    }
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct Sphere<Center, Radius>
where
    Center: for<'a> BinRead<Args<'a> = ()>, // Use an explicit lifetime 'a
    Radius: for<'a> BinRead<Args<'a> = ()>, // Use an explicit lifetime 'a
{
    pub center: Center,
    pub radius: Radius,
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct MinMax<T>
where
    T: for<'a> BinRead<Args<'a> = ()>, // Use an explicit lifetime 'a
{
    min: T,
    max: T,
}

#[binread]
#[derive(Debug)]
#[br(little)]
#[derive(Default)]
pub struct EdgeIndicesShort {
    #[br(dbg)]
    pub westVertexCount: u32,
    #[br(count = westVertexCount)]
    pub westIndices: Vec<u16>,

    #[br(dbg)]
    pub southVertexCount: u32,
    #[br(count = southVertexCount)]
    pub southIndices: Vec<u16>,

    #[br(dbg)]
    pub eastVertexCount: u32,
    #[br(count = eastVertexCount)]
    pub eastIndices: Vec<u16>,

    #[br(dbg)]
    pub northVertexCount: u32,
    #[br(count = northVertexCount)]
    pub northIndices: Vec<u16>,
}

#[binread]
#[derive(Debug)]
#[br(little)]
#[derive(Default)]
pub struct EdgeIndicesLong {
    pub westVertexCount: u32,
    #[br(count = westVertexCount)]
    pub westIndices: Vec<u32>,

    pub southVertexCount: u32,
    #[br(count = southVertexCount)]
    pub southIndices: Vec<u32>,

    pub eastVertexCount: u32,
    #[br(count = eastVertexCount)]
    pub eastIndices: Vec<u32>,

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
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i))]
    pub u: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i))]
    pub v: Vec<i32>,

    #[br(count = vertex_count)]
    #[br(map = |i:  Vec<u16>| vec_zigzag_decode(i))]
    pub height: Vec<i32>,

    #[br(dbg)]
    pub triangleCount: u32,

    #[br(dbg)]
    #[br(align_before = 4)]
    #[br(count = triangleCount * 3)]
    #[br(if(vertex_count > 65536) )]
    #[br(map = |i:  Vec<u32>| i.high_watermark_decode())]
    pub index_data_long: Vec<u32>,

    #[br(align_before = 2)]
    #[br(count = triangleCount * 3)]
    #[br(if(vertex_count <= 65536))]
    #[br(map = |i:  Vec<u16>| i.high_watermark_decode())]
    pub index_data_short: Vec<u16>,

    #[br(align_before = 4)]
    #[br(if(vertex_count > 65536) )]
    pub edge_indices_long: EdgeIndicesLong,

    #[br(align_before = 2)]
    #[br(if(vertex_count <= 65536) )]
    pub edge_indices_short: EdgeIndicesShort,
}

trait HighWatermarkDecode<'a, T> {
    fn high_watermark_decode(&'a self) -> Vec<T>;
}

impl<'a, T> HighWatermarkDecode<'a, T> for Vec<T>
where
    T: 'a + Copy + Sub<Output = T> + From<u8> + PartialEq + Add<Output = T>,
{
    fn high_watermark_decode(&'a self) -> Vec<T> {
        let mut decoded: Vec<T> = Vec::with_capacity(self.len());
        let mut highest: T = T::from(0u8);

        for &i in self.iter() {
            decoded.push(highest - i);
            if i == T::from(0u8) {
                highest = highest + T::from(1u8);
            }
        }
        decoded
        //decoded.iter().tuple_windows::<Triangle<T>>()
    }
}

impl VertexData {
    pub fn as_uvh(&self) -> Vec<Vertex3<&i32>> {
        izip!(&self.u, &self.v, &self.height).collect()
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

pub type Center = CartesianPoint3<f64>;
pub type BoundingSphereCenter = Sphere<CartesianPoint3<f64>, f64>;
pub type HorizonOcclusionPoint = CartesianPoint3<f64>;
pub type Height = MinMax<f32>;

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct QuantizedMeshHeader {
    pub center: Center,
    pub height: Height,
    pub bounding_sphere_center: BoundingSphereCenter,
    pub horizon_occlusion_point: HorizonOcclusionPoint,
}

#[binread]
#[derive(Debug)]
#[br(little)]

pub struct QuantizedMesh {
    pub header: QuantizedMeshHeader,
    pub vertex_data: VertexData,
}

pub fn localecef2geodetic(
    pt: CartesianPoint3<f64>,
    r_ellips: Ellipsoid,
) -> GeodeticPoint3<f64, f64> {
    let (major, minor, _, _) = r_ellips.parameters();

    let r = (pt.x * pt.x + pt.y * pt.y + pt.z * pt.z).sqrt();
    let e = (major * major - minor * minor).sqrt();
    let var = r * r - e * e;
    let u = (0.5 * var + 0.5 * (var * var + 4.0 * e * e * pt.z * pt.z).sqrt()).sqrt();

    let q = (pt.x * pt.x + pt.y * pt.y).sqrt();
    let hu_e = (u * u + e * e).sqrt();
    let mut beta = (hu_e / u * pt.z / q).atan();

    let eps = ((minor * u - major * hu_e + e * e) * beta.sin())
        / (major * hu_e / beta.cos() - e * e * beta.cos());
    beta += eps;

    let lat: Radian<f64> = (major / minor * beta.tan()).atan();
    let lon: Radian<f64> = pt.y.atan2(pt.x);

    let v1 = pt.z - minor * beta.sin();
    let v2 = q - major * beta.cos();
    let alt: Meter<f64>;

    let inside = (pt.x * pt.x / major / major)
        + (pt.y * pt.y / major / major)
        + (pt.z * pt.z / minor / minor)
        < 1.0;
    if inside {
        alt = -(v1 * v1 + v2 * v2).sqrt();
    } else {
        alt = (v1 * v1 + v2 * v2).sqrt();
    };

    GeodeticPoint3 { lat, lon, alt }
}
