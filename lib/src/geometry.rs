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

        let a = T::from_f64(ellipsoid.semi_major_axis).unwrap();
        let b = T::from_f64(ellipsoid.semi_minor_axis).unwrap();
        let e2 = (a * a - b * b) / (a * a);
        let e_prime2 = (a * a - b * b) / (b * b);

        let p = (x * x + y * y).sqrt();
        let theta = (z * a).atan2(p * b);

        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        let lat = ((z + e_prime2 * b * sin_theta * sin_theta * sin_theta) / (p - e2 * a * cos_theta * cos_theta * cos_theta)).atan();
        let lon = y.atan2(x);

        let sin_lat = lat.sin();
        let n = a / (T::one() - e2 * sin_lat * sin_lat).sqrt();
        let alt = p / lat.cos() - n;


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
        &self.0.x
    }

    pub fn lon(&self) -> &T {
        &self.0.y
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

        let n = a / ((T::one() -  T::from_f64(ellipsoid.eccentricity_squared).unwrap() * lat.sin().powi(2)).sqrt()); // Radius of curvature in the prime vertical

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
        geodetic_surface_normal(self, p)
    }
}

// Shared Functions
pub fn geodetic_surface_normal(ellipsoid: &Ellipsoid, p: &CartesianPoint3<f64>) -> Vector3<f64> {
    let normal = Vector3::new(
        p.0.x / (ellipsoid.semi_major_axis * ellipsoid.semi_major_axis),
        p.0.y / (ellipsoid.semi_major_axis * ellipsoid.semi_major_axis),
        p.0.z / (ellipsoid.semi_minor_axis * (ellipsoid.semi_minor_axis)),
    );
    // Normalize the vector to get the unit normal
    normal.normalize()
}
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
       0.0,    0.0,    0.0, 1.0,
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

 #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn test_cartesian_point3_new() {
            let point = CartesianPoint3::new(1.0, 2.0, 3.0);
            assert_eq!(point.0, Point3::new(1.0, 2.0, 3.0));
        }

        #[test]
        fn test_geodetic_point3_new() {
            let point = GeodeticPoint3::new(1.0, 2.0, 3.0);
            assert_eq!(point.0, Point3::new(1.0, 2.0, 3.0));
        }

        #[test]
        fn test_geodetic_point3_from_degrees() {
            let point = GeodeticPoint3::from_degrees(180.0, 90.0, 100.0);
            assert_eq!(point.0, Point3::new(std::f64::consts::PI, std::f64::consts::PI / 2.0, 100.0));
        }

        #[test]
        fn test_geodetic_point3_to_degrees() {
            let point = GeodeticPoint3::new(std::f64::consts::PI, std::f64::consts::PI / 2.0, 100.0);
            let degrees_point = point.to_degrees();
            assert_eq!(degrees_point.0, Point3::new(180.0, 90.0, 100.0));
        }

        #[test]
        fn test_geodetic_point3_to_ecef() {
            let ellipsoid = Ellipsoid::wgs84();
            let geodetic_point = GeodeticPoint3::new(0.0, 0.0, 0.0);
            let ecef_point = geodetic_point.to_ecef(&ellipsoid);
            assert_eq!(ecef_point.0, Point3::new(ellipsoid.semi_major_axis, 0.0, 0.0));
        }

        #[test]
        fn test_rectangle_height() {
            let rect = Rectangle {
                lower_left: Point2::new(0.0, 0.0),
                upper_right: Point2::new(1.0, 1.0),
            };
            assert_eq!(rect.height(), 1.0);
        }

        #[test]
        fn test_rectangle_width() {
            let rect = Rectangle {
                lower_left: Point2::new(0.0, 0.0),
                upper_right: Point2::new(1.0, 1.0),
            };
            assert_eq!(rect.width(), 1.0);
        }

        #[test]
        fn test_geodetic_rectangle_to_degrees() {
            let rect = GeodeticRectangle {
                lower_left: GeodeticPoint2::from_degrees(0.0, 0.0),
                upper_right: GeodeticPoint2::from_degrees(1.0, 1.0),
            };
            let degrees_rect = rect.to_degrees();
            assert_eq!(degrees_rect, ((0.0, 0.0), (1.0, 1.0)));
        }

        #[test]
        fn test_triangle_bounding_rect() {
            let triangle = Triangle {
                vertices: [
                    Point2::new(0.0, 0.0),
                    Point2::new(1.0, 1.0),
                    Point2::new(0.5, 0.5),
                ],
            };
            let bounding_rect = triangle.bounding_rect();
            assert_eq!(bounding_rect.lower_left, Point2::new(0.0, 0.0));
            assert_eq!(bounding_rect.upper_right, Point2::new(1.0, 1.0));
        }

        #[test]
        fn test_calculate_enu_to_ecef_rotation_matrix() {
            let ellipsoid = Ellipsoid::wgs84();
            let ecef_position = CartesianPoint3::new(ellipsoid.semi_major_axis, 0.0, 0.0);
            let rotation_matrix = calculate_enu_to_ecef_rotation_matrix(&ecef_position, &ellipsoid);
            println!("{:?}", rotation_matrix);
            assert_eq!(rotation_matrix,  Matrix4::new(
                0.0, 0.0, 1.0, 0.0,
                1.0, 0.0, 0.0, 0.0,
                0.0, 1.0, 0.0, 0.0,
                6378137.0, 0.0, 0.0, 1.0
            ));
        }

        #[test]
        fn test_lerp() {
            let result = lerp(&0.0, &10.0, &0.5);
            assert_eq!(result,   5.0);
        }

        #[test]
        fn test_radians_to_degrees() {
            let radians = std::f64::consts::PI;
            let degrees = radians_to_degrees(radians);
            assert_eq!(degrees, 180.0);
        }

        #[test]
        fn test_degrees_to_radians() {
            let degrees = 180.0;
            let radians = degrees_to_radians(degrees);
            assert_eq!(radians, std::f64::consts::PI);
        }
    }