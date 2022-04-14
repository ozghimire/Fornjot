use std::cmp::max;

use fj_math::{Point, Scalar, Transform, Vector};

use crate::algorithms::Tolerance;

/// An arc
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
pub struct Arc {
    /// The center point of the circle the arc is on
    pub center: Point<3>,

    /// The radius of the circle the arc is on
    ///
    /// The radius is represented by a vector that points from the center to the
    /// circumference. The point on the circumference that it points to defines
    /// the origin of the arc's 1-dimensional curve coordinate system.
    pub radius: Vector<2>,

    /// The length of the arc in radians
    ///
    /// If the length is positive, the arc is oriented counter-clockwise. If it
    /// is negative, the arc is oriented clockwise.
    pub length: Scalar,
}

impl Arc {
    /// Access the origin of the curve's coordinate system
    pub fn origin(&self) -> Point<3> {
        self.center
    }

    /// Create a new instance that is transformed by `transform`
    #[must_use]
    pub fn transform(self, transform: &Transform) -> Self {
        let radius = self.radius.to_xyz();
        let radius = transform.transform_vector(&radius);
        let radius = radius.xy();

        Self {
            center: transform.transform_point(&self.center),
            radius,
            length: self.length,
        }
    }

    /// Convert a point in model coordinates to curve coordinates
    ///
    /// Converts the provided point into curve coordinates between `0.`
    /// (inclusive) and `PI * 2.` (exclusive).
    ///
    /// Projects the point onto the circle before computing curve coordinate,
    /// ignoring the radius. This is done to make this method robust against
    /// floating point accuracy issues.
    ///
    /// Callers are advised to be careful about the points they pass, as the
    /// point not being on the curve, intentional or not, will not result in an
    /// error.
    pub fn point_model_to_curve(&self, point: &Point<3>) -> Point<1> {
        let v = point - self.center;
        let atan = Scalar::atan2(v.y, v.x);
        let coord = if atan >= Scalar::ZERO {
            atan
        } else {
            atan + Scalar::PI * 2.
        };
        Point::from([coord])
    }

    /// Convert a point on the curve into model coordinates
    pub fn point_curve_to_model(&self, point: &Point<1>) -> Point<3> {
        self.center + self.vector_curve_to_model(&point.coords)
    }

    /// Convert a vector on the curve into model coordinates
    pub fn vector_curve_to_model(&self, vector: &Vector<1>) -> Vector<3> {
        let radius = self.radius.magnitude();
        let angle = vector.t;

        let (sin, cos) = angle.sin_cos();

        let x = cos * radius;
        let y = sin * radius;

        Vector::from([x, y, Scalar::ZERO])
    }

    /// Approximate the arc
    ///
    /// `tolerance` specifies how much the approximation is allowed to deviate
    /// from the arc.
    pub fn approx(&self, tolerance: Tolerance, out: &mut Vec<Point<3>>) {
        let radius = self.radius.magnitude();

        // To approximate the circle, we use a regular polygon for which
        // the circle is the circumscribed circle. The `tolerance`
        // parameter is the maximum allowed distance between the polygon
        // and the circle. This is the same as the difference between
        // the circumscribed circle and the incircle.

        let n = Self::number_of_vertices(tolerance, radius, self.length);

        for i in 0..n {
            let angle = self.length / n as f64 * i as f64;
            let point = self.point_curve_to_model(&Point::from([angle]));
            out.push(point);
        }
    }

    fn number_of_vertices(
        tolerance: Tolerance,
        radius: Scalar,
        length: Scalar,
    ) -> u64 {
        let length = length.abs();

        let n =
            (length / 2. / (Scalar::ONE - (tolerance.inner() / radius)).acos())
                .ceil()
                .into_u64();

        if length >= Scalar::PI * 2. {
            max(n, 3)
        } else {
            max(n, 2)
        }
    }
}

#[cfg(test)]
mod tests {
    use std::f64::consts::{FRAC_PI_2, PI};

    use fj_math::{Point, Scalar, Vector};

    use crate::algorithms::Tolerance;

    use super::Arc;

    #[test]
    fn point_model_to_curve() {
        let circle = Arc {
            center: Point::from([1., 2., 3.]),
            radius: Vector::from([1., 0.]),
            length: Scalar::PI * 2.,
        };

        assert_eq!(
            circle.point_model_to_curve(&Point::from([2., 2., 3.])),
            Point::from([0.]),
        );
        assert_eq!(
            circle.point_model_to_curve(&Point::from([1., 3., 3.])),
            Point::from([FRAC_PI_2]),
        );
        assert_eq!(
            circle.point_model_to_curve(&Point::from([0., 2., 3.])),
            Point::from([PI]),
        );
        assert_eq!(
            circle.point_model_to_curve(&Point::from([1., 1., 3.])),
            Point::from([FRAC_PI_2 * 3.]),
        );
    }

    #[test]
    fn number_of_vertices() {
        // full circle
        verify_result(50., 100., Scalar::PI * 2., 3);
        verify_result(10., 100., Scalar::PI * 2., 7);
        verify_result(1., 100., Scalar::PI * 2., 23);

        // half circle
        verify_result(50., 100., Scalar::PI, 2);
        verify_result(10., 100., Scalar::PI, 4);
        verify_result(1., 100., Scalar::PI, 12);

        fn verify_result(
            tolerance: impl Into<Tolerance>,
            radius: impl Into<Scalar>,
            length: Scalar,
            n: u64,
        ) {
            let tolerance = tolerance.into();
            let radius = radius.into();

            assert_eq!(n, Arc::number_of_vertices(tolerance, radius, length));

            assert!(calculate_error(radius, length, n) <= tolerance.inner());
            if n > 3 {
                assert!(
                    calculate_error(radius, length, n - 1) >= tolerance.inner()
                );
            }
        }

        fn calculate_error(radius: Scalar, length: Scalar, n: u64) -> Scalar {
            radius - radius * (length / 2. / Scalar::from_u64(n)).cos()
        }
    }
}
