use std::{collections::HashMap, f64::consts::PI};

use nalgebra::vector;
use parry3d_f64::math::Isometry;

use crate::{
    debug::DebugInfo,
    kernel::{
        algorithms::{approximation::Approximation, transform::transform_face},
        geometry::{surfaces::Swept, Curve, Line, Surface},
        topology::{
            edges::{Cycle, Edge, Edges},
            faces::{Face, Faces},
            Shape,
        },
    },
    math::{Aabb, Scalar, Transform, Vector},
};

use super::ToShape;

impl ToShape for fj::Sweep {
    fn to_shape(&self, tolerance: Scalar, debug_info: &mut DebugInfo) -> Shape {
        let mut shape = Shape::new();

        let original_shape = self.shape.to_shape(tolerance, debug_info);

        let path = Vector::from([0., 0., self.length]);
        let rotation = Isometry::rotation(vector![PI, 0., 0.]).into();
        let translation = Isometry::translation(0.0, 0.0, self.length).into();

        let mut bottom_faces = Vec::new();
        let mut top_faces = Vec::new();
        let mut side_faces = Vec::new();

        for face in original_shape.faces.0 {
            // This only works for faces that are symmetric to the x-axis.
            //
            // See issue:
            // https://github.com/hannobraun/Fornjot/issues/230
            bottom_faces.push(transform_face(&face, &rotation, &mut shape));

            top_faces.push(transform_face(&face, &translation, &mut shape));
        }

        for cycle in original_shape.edges.cycles {
            if cycle.edges.len() == 1 {
                // If there's only one edge, it must be a continuous edge that
                // connects to itself. By sweeping that, we create a continuous
                // face.
                //
                // Continuous faces aren't currently supported by the
                // approximation code, and hence can't be triangulated. To
                // address that, we fall back to the old and almost obsolete
                // triangle representation to create the face.
                //
                // This is the last piece of code that still uses the triangle
                // representation.

                let approx = Approximation::for_cycle(&cycle, tolerance);

                let mut quads = Vec::new();
                for segment in approx.segments {
                    let [v0, v1] = segment.points();
                    let [v3, v2] = {
                        let segment =
                            Transform::translation(0., 0., self.length)
                                .transform_segment(&segment);
                        segment.points()
                    };

                    quads.push([v0, v1, v2, v3]);
                }

                let mut side_face = Vec::new();
                for [v0, v1, v2, v3] in quads {
                    side_face.push([v0, v1, v2].into());
                    side_face.push([v0, v2, v3].into());
                }

                side_faces.push(Face::Triangles(side_face));
            } else {
                // If there's no continuous edge, we can create the non-
                // continuous faces using boundary representation.

                let mut side_edges = HashMap::new();

                for bottom_edge in cycle.edges {
                    let top_edge = bottom_edge.clone().transform(&translation);

                    let surface = Surface::Swept(Swept {
                        curve: bottom_edge.curve.clone(),
                        path,
                    });

                    let side_edges = bottom_edge.vertices.map(|vertices| {
                        vertices.map(|vertex| {
                            let edge =
                                side_edges.entry(vertex).or_insert_with(|| {
                                    let line = Line {
                                        origin: vertex.point().canonical(),
                                        direction: path,
                                    };

                                    let a = vertex;
                                    // TASK: This is wrong. The vertex created
                                    //       here already exists in `top_edges`.
                                    let b = shape
                                        .vertices()
                                        .create(line.origin + line.direction);

                                    let curve = Curve::Line(line);

                                    let vertices = Some([a.to_canonical(), b]);

                                    Edge::new(curve, vertices)
                                });

                            edge.clone()
                        })
                    });

                    let face = match side_edges {
                        Some([a, b]) => Face::Face {
                            surface,
                            edges: Edges::single_cycle([
                                bottom_edge,
                                top_edge,
                                a,
                                b,
                            ]),
                        },
                        None => Face::Face {
                            surface,
                            edges: Edges {
                                cycles: vec![
                                    Cycle {
                                        edges: vec![bottom_edge],
                                    },
                                    Cycle {
                                        edges: vec![top_edge],
                                    },
                                ],
                            },
                        },
                    };

                    side_faces.push(face);
                }
            }
        }

        let mut faces = Vec::new();
        faces.extend(bottom_faces);
        faces.extend(top_faces);
        faces.extend(side_faces);

        shape.faces = Faces(faces);

        shape
    }

    fn bounding_volume(&self) -> Aabb<3> {
        let mut aabb = self.shape.bounding_volume();
        aabb.max.z = self.length.into();
        aabb
    }
}
