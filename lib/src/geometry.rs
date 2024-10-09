use nalgebra::{Matrix4, Point2, Point3, RealField, Scalar, SimdPartialOrd, Vector3};
use std::clone::Clone;
use std::cmp::PartialEq;
use std::fmt::Debug;
use std::ops::Sub;

// nalgebra Wrappers
#[derive(Debug)]
pub struct CartesianPoint3<T>(pub Point3<T>)
where
    T: Copy + RealField;

impl<T> CartesianPoint3<T>
where
    T: Copy + RealField ,
{
    // Constructor for ECEFPoint3
    pub fn new(x: T, y: T, z: T) -> Self {
        CartesianPoint3(Point3::new(x, y, z))
    }

    pub fn to_geodetic(&self, ellipsoid: &Ellipsoid) -> GeodeticPoint3<T> {
        let x = self.0.x;
        let y = self.0.y;
        let z = self.0.z;

        let major = T::from_f64(ellipsoid.semi_major_axis).unwrap();
        let minor = T::from_f64(ellipsoid.semi_minor_axis).unwrap();

        let r = (x * x + y * y + z * z).sqrt();
        let e = (major * major - minor * minor).sqrt();
        let var = r * r - e * e;
        let u = (T::from_f64(0.5).unwrap() * var
            + T::from_f64(0.5).unwrap()
                * (var * var + T::from_f64(4.0).unwrap() * e * e * z * z).sqrt())
        .sqrt();

        let q = (x * x + y * y).sqrt();
        let hu_e = (u * u + e * e).sqrt();
        let mut beta = (hu_e / u * z / q).atan();

        let eps = ((minor * u - major * hu_e + e * e) * beta.sin())
            / (major * hu_e / beta.cos() - e * e * beta.cos());
        beta += eps;

        let lat = (major / minor * beta.tan()).atan();
        let lon = y.atan2(x);

        let v1 = z - minor * beta.sin();
        let v2 = q - major * beta.cos();
        let inside =
            (x * x / major / major) + (y * y / major / major) + (z * z / minor / minor)
                < T::one();
        let alt = if inside {
            -(v1 * v1 + v2 * v2).sqrt()
        } else {
            (v1 * v1 + v2 * v2).sqrt()
        };

        GeodeticPoint3::new(lat, lon, alt) // Direct construction
    }
}

impl<T> Sub for CartesianPoint3<T>
where
    T: Copy + RealField + Sub<Output = T>,
{
    type Output = CartesianPoint3<T>;

    fn sub(self, other: CartesianPoint3<T>) -> Self::Output {
        CartesianPoint3((self.0 - other.0).into()) // Use nalgebra's subtraction
    }
}

//  GeodeticPoint3
#[derive(Debug)]
pub struct GeodeticPoint3<T>(pub Point3<T>)
where
    T: Copy + Scalar;

impl<T> GeodeticPoint3<T>
where
    T: Copy + RealField,
{
    // Constructor for GeodeticPoint3
    pub fn new(lat: T, lon: T, alt: T) -> Self {
        GeodeticPoint3(Point3::new(lat, lon, alt))
    }

    pub fn lat(&self) -> &T {
        &self.0.y
    }

    pub fn lon(&self) -> &T {
        &self.0.x
    }

    pub fn height(&self) -> &T {
        &self.0.z
    }

    pub fn from_degrees(lat: T, lon: T, alt: T) -> Self {
        GeodeticPoint3(Point3::new(
            degrees_to_radians(lat),
            degrees_to_radians(lon),
            alt,
        ))
    }
    pub fn to_degrees(&self) -> GeodeticPoint3<T> {
        GeodeticPoint3::new(
            radians_to_degrees(*self.lat()),
            radians_to_degrees(*self.lon()),
            *self.height(),
        )
    }

    pub fn to_ecef(&self, ellipsoid: &Ellipsoid) -> CartesianPoint3<T> {
        let lat = self.lat();
        let lon = self.lon();
        let height = self.height();

        let a = T::from_f64(ellipsoid.semi_major_axis).unwrap();
        let b = T::from_f64(ellipsoid.semi_minor_axis).unwrap();
        let e2 = (a * a - b * b) / (a * a); // Square of eccentricity

        let n = a / ((T::one() - e2 * lat.sin().powi(2)).sqrt()); // Radius of curvature in the prime vertical

        let x = (n + *height) * lat.cos() * lon.cos();
        let y = (n + *height) * lat.cos() * lon.sin();
        let z = (n * (T::one() - e2) + *height) * lat.sin();

        CartesianPoint3(Point3::new(x, y, z))
    }
}


// GeodeticPoint2
#[derive(Debug, Default)]
pub struct GeodeticPoint2<T>(pub Point2<T>)
where
    T: Copy + RealField;

impl<T> GeodeticPoint2<T>
where
    T: Copy + RealField,
{
    pub fn new(lat: T, lon: T) -> Self {
        GeodeticPoint2(Point2::new(lat, lon))
    }

    pub fn from_degrees(lat: T, lon: T) -> Self {
        GeodeticPoint2::new(degrees_to_radians(lat), degrees_to_radians(lon))
    }
    pub fn to_degrees(&self) -> (T,T) {
        (
            radians_to_degrees(*self.lat()),
            radians_to_degrees(*self.lon()),
        )
    }
    pub fn lat(&self) -> &T {
        &self.0.y
    }

    pub fn lon(&self) -> &T {
        &self.0.x
    }
}



// Polygons

// Rectangle
#[derive(Debug)]
pub struct Rectangle<T>
where
    T: Clone + Debug + PartialEq + 'static,
{
    pub lower_left: Point2<T>,
    pub upper_right: Point2<T>,
}

impl<T> Rectangle<T>
where
    T: RealField + Copy + Clone + Debug + PartialEq + Sub<Output = T>,
{
    pub fn height(&self) -> T {
        self.upper_right.y - self.lower_left.y
    }

    pub fn width(&self) -> T {
        self.upper_right.x - self.lower_left.x
    }
}



#[derive(Debug, Default)]
pub struct GeodeticRectangle<T>
where
    T: Copy + RealField + Clone + Debug + PartialEq + 'static,
{
    pub lower_left: GeodeticPoint2<T>,
    pub upper_right: GeodeticPoint2<T>,
}

impl<T> GeodeticRectangle<T>
where
    T: RealField + Copy + Clone + Debug + PartialEq + Sub<Output = T>,
{
    pub fn height(&self) -> T {
        self.upper_right.0.y - self.lower_left.0.y
    }

    pub fn width(&self) -> T {
        self.upper_right.0.x - self.lower_left.0.x
    }

     /// Converts the rectangle's corners from radians to degrees.
     pub fn to_degrees(&self) -> ((T, T), (T, T)) {
        (
            self.lower_left.to_degrees(),
            self.upper_right.to_degrees(),
        )
    }
}

// Triangle

pub struct Triangle<T>
where
    T: Clone + Debug + PartialEq + 'static,
{
    pub vertices: [Point2<T>; 3], // Store three vertices
}

impl<T> Triangle<T>
where
    T: Copy + Scalar + Clone + Debug + SimdPartialOrd,
{
    // Ensure Point supports required operations{
    pub fn bounding_rect(&self) -> Rectangle<T> {
        // componentwise evaluation of the triangle bounding box
        let (mut lower_left, mut upper_right) = self.vertices[0].inf_sup(&self.vertices[1]);
        lower_left = lower_left.inf(&self.vertices[2]);
        upper_right = upper_right.sup(&self.vertices[2]);

        Rectangle {
            lower_left,
            upper_right,
        }
    }
}

/// Ellipsoids
#[derive(Debug, Default)]
pub struct Ellipsoid {
    pub semi_major_axis: f64, // Semi-major axis (equatorial radius) in meters
    pub semi_minor_axis: f64, // Semi-major axis (equatorial radius) in meters
    pub flattening: f64,      // Semi-minor axis (polar radius) in metsers
    pub eccentricity_squared: f64,
}

impl Ellipsoid {
    /// Create a WGS-84 Ellipsoid
    pub fn wgs84() -> Self {
        let semi_major_axis = 6378137.0;
        let flattening = 1.0 / 298.257223563;
        let semi_minor_axis = semi_major_axis * (1.0 - flattening);
        let eccentricity_squared = flattening * (2.0 - flattening);

        Self {
            semi_major_axis,
            semi_minor_axis,
            flattening,
            eccentricity_squared,
        }
    }

    /// Calculate the surface normal vector at the given ECEF position
    pub fn geodetic_surface_normal(&self, p: &CartesianPoint3<f64>) -> Vector3<f64> {
        let normal = Vector3::new(
            p.0.x / (self.semi_major_axis * self.semi_major_axis),
            p.0.y / (self.semi_major_axis * self.semi_major_axis),
            p.0.z / (self.semi_minor_axis * (self.semi_minor_axis)),
        );
        // Normalize the vector to get the unit normal
        normal.normalize()
    }
}

// Shared Functions

/// Calculate the 4x4 ENU to ECEF rotation matrix for a given ECEF reference point
pub fn calculate_enu_to_ecef_rotation_matrix(
    ecef_position: &CartesianPoint3<f64>,
    ellipsoid: &Ellipsoid,
) -> Matrix4<f64> {
    let up = ellipsoid.geodetic_surface_normal(ecef_position);
    let mut east: Vector3<f64> = Vector3::new(ecef_position.0.y, ecef_position.0.x, 0.0);
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
        ecef_position.0.x,
        ecef_position.0.y,
        ecef_position.0.z,
        1.0,
    )
}

// Linear interpolation
pub fn lerp(min_value: &f64, max_value: &f64, t: &f64) -> f64 {
    min_value + t * (max_value - min_value)
}

// Generic radians/degrees transformation

fn radians_to_degrees<T>(radians: T) -> T
where
    T: RealField, // Ensures T supports real number operations
{
    let pi = T::pi(); // Get the value of pi for the generic type T
    let degrees_per_radian = T::from_f64(180.0).unwrap() / pi; // 180 / pi

    radians * degrees_per_radian // Convert radians to degrees
}

fn degrees_to_radians<T>(degrees: T) -> T
where
    T: RealField, // Ensures T supports real number operations
{
    let pi = T::pi(); // Get the value of pi for the generic type T
    let degrees_per_radian = T::from_f64(180.0).unwrap() / pi; // 180 / pi

    degrees / degrees_per_radian // Convert radians to degrees
}
