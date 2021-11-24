use std::f32::consts::FRAC_PI_4;

use nalgebra::{Isometry3, Perspective3, Rotation, Translation};

#[derive(Debug)]
pub struct Transform {
    pub rotation: Rotation<f32, 3>,
    pub translation: Translation<f32, 2>,
    pub distance: f32,
}

impl Transform {
    pub fn new(initial_distance: f32) -> Self {
        Self {
            rotation: Rotation::identity(),
            translation: Translation::identity(),
            distance: initial_distance,
        }
    }

    pub fn to_native(&self, aspect_ratio: f32) -> NativeTransform {
        let projection = Perspective3::new(
            aspect_ratio,
            FIELD_OF_VIEW,
            NEAR_PLANE,
            FAR_PLANE,
        );

        let transform = projection.to_projective() * self.view_transform();

        let mut native = [0.0; 16];
        native.copy_from_slice(transform.matrix().data.as_slice());

        native
    }

    pub fn to_normals_transform(&self) -> NativeTransform {
        let transform =
            self.view_transform().inverse().to_homogeneous().transpose();

        let mut native = [0.0; 16];
        native.copy_from_slice(transform.data.as_slice());

        native
    }

    fn view_transform(&self) -> Isometry3<f32> {
        Isometry3::from_parts(
            Translation::from([
                self.translation.x,
                self.translation.y,
                -self.distance,
            ]),
            self.rotation.into(),
        )
    }
}

pub type NativeTransform = [f32; 16];

pub const NEAR_PLANE: f32 = 0.1;
pub const FAR_PLANE: f32 = 1000.0;
pub const FIELD_OF_VIEW: f32 = FRAC_PI_4; // 45 degrees
