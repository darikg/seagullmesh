#include "seagullmesh.hpp"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/angle_and_area_smoothing.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/tangential_relaxation.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/Adaptive_sizing_field.h>
#include <CGAL/Polygon_mesh_processing/Uniform_sizing_field.h>

typedef std::vector<V>                                      Verts;
typedef std::vector<F>                                      Faces;
typedef Mesh3::Property_map<V, Point3>                      VertPoint;
typedef Mesh3::Property_map<V, double>                      VertDouble;
typedef Mesh3::Property_map<V, bool>                        VertBool;
typedef Mesh3::Property_map<E, bool>                        EdgeBool;
typedef Mesh3::Property_map<F, F>                           FaceMap;

namespace PMP = CGAL::Polygon_mesh_processing;

typedef PMP::Principal_curvatures_and_directions<Kernel>    PrincipalCurvDir;
typedef Mesh3::Property_map<V, PrincipalCurvDir>            VertPrincipalCurvDir;
typedef PMP::Adaptive_sizing_field<Mesh3, VertPoint>        AdaptiveSizingField;
typedef PMP::Uniform_sizing_field<Mesh3, VertPoint>         UniformSizingField;


struct VertexPointMapWrapper {
    // Used for tracking which verts get moved during remesh, etc
    using key_type = V;
    using value_type = Point3;
    using reference = Point3&;
    using category = boost::read_write_property_map_tag;

    VertPoint& points;
    VertBool& touched;

    VertexPointMapWrapper(VertPoint& p, VertBool& t) : points(p), touched(t) {}

    friend Point3& get (const VertexPointMapWrapper& map, V v) { return map.points[v]; }
    friend void put (const VertexPointMapWrapper& map, V v, const Point3& point) {
        map.points[v] = point;
        map.touched[v] = true;
    }
};

template <typename Vpm, typename SizingField>
void define_isotropic_remeshing(py::module &m) {
    m.def("remesh", [](
        Mesh3& mesh,
        const Faces& faces,
        SizingField& sizing_field,
        unsigned int n_iter,
        bool protect_constraints,
        VertBool& vertex_is_constrained_map,
        EdgeBool& edge_is_constrained_map,
        Vpm& vertex_point_map
    ) {
        auto params = PMP::parameters::
            number_of_iterations(n_iter)
            .vertex_point_map(vertex_point_map)
            .protect_constraints(protect_constraints)
            .vertex_is_constrained_map(vertex_is_constrained_map)
            .edge_is_constrained_map(edge_is_constrained_map)
        ;
        PMP::isotropic_remeshing(faces, sizing_field, mesh, params);
    });
}

void init_meshing(py::module &m) {

    py::module sub = m.def_submodule("meshing");
    define_isotropic_remeshing<VertPoint,               UniformSizingField> (sub);
    define_isotropic_remeshing<VertPoint,               AdaptiveSizingField>(sub);
    define_isotropic_remeshing<VertexPointMapWrapper,   UniformSizingField> (sub);
    define_isotropic_remeshing<VertexPointMapWrapper,   UniformSizingField>(sub);

    sub.def("fair", [](Mesh3& mesh, const Verts& verts, const unsigned int fairing_continuity) {
            // A value controling the tangential continuity of the output surface patch.
            // The possible values are 0, 1 and 2, refering to the C0, C1 and C2 continuity.
            auto params = PMP::parameters::fairing_continuity(fairing_continuity);
            bool success = PMP::fair(mesh, verts, params);
            if (!success) {
                throw std::runtime_error("Fairing failed");
            }
        })
        .def("refine", [](Mesh3& mesh, const Faces& faces, double density) {
            std::vector<V> new_verts;
            std::vector<F> new_faces;
            auto params = PMP::parameters::density_control_factor(density);
            PMP::refine(mesh, faces, std::back_inserter(new_faces), std::back_inserter(new_verts), params);
            return std::make_tuple(new_verts, new_faces);
        })
        .def("smooth_angle_and_area", [](
            Mesh3& mesh, 
            const std::vector<F>& faces, 
            unsigned int n_iter,
            bool use_area_smoothing,
            bool use_angle_smoothing,
            bool use_safety_constraints,
            bool do_project, 
            VertBool& vertex_is_constrained_map,
            EdgeBool& edge_is_constrained_map
        ) {
            auto params = PMP::parameters::
                number_of_iterations(n_iter)
                .use_area_smoothing(use_area_smoothing)
                .use_angle_smoothing(use_angle_smoothing)
                .use_safety_constraints(use_safety_constraints)
                .do_project(do_project)
                .vertex_is_constrained_map(vertex_is_constrained_map)
                .edge_is_constrained_map(edge_is_constrained_map)
            ;
            PMP::angle_and_area_smoothing(faces, mesh, params);
        })
        .def("tangential_relaxation", [](
            Mesh3& mesh,
            const std::vector<V>& verts,
            unsigned int n_iter,
            bool relax_constraints,
            VertBool& vertex_is_constrained_map,
            EdgeBool& edge_is_constrained_map
        ) {
            auto params = PMP::parameters::
                number_of_iterations(n_iter)
                .relax_constraints(relax_constraints)
                .vertex_is_constrained_map(vertex_is_constrained_map)
                .edge_is_constrained_map(edge_is_constrained_map)
            ;
            PMP::tangential_relaxation(verts, mesh, params);
        })
        .def("smooth_shape", [](
            Mesh3& mesh, 
            const std::vector<F>& faces, 
            const double time, 
            unsigned int n_iter,
            VertBool& vertex_is_constrained_map
        ) {
            auto params = PMP::parameters::
                number_of_iterations(n_iter)
                .vertex_is_constrained_map(vertex_is_constrained_map)
            ;
            PMP::smooth_shape(faces, mesh, time, params);
        })
        .def("does_self_intersect", [](Mesh3& mesh) {
            return PMP::does_self_intersect(mesh);
        })
        .def("remesh_planar_patches", [](
                const Mesh3& mesh,
                EdgeBool& edge_is_constrained_map,
                // FaceMap& face_patch_map,
                float cosine_of_maximum_angle
            ) {
            // TODO parameters.face_patch_map wants propertymap<F, std::size_t>
            // But we've only exposed pmap<F, int>
            // Will doing both signed and unsigned ints be too complicated?
            auto params = PMP::parameters::
                edge_is_constrained_map(edge_is_constrained_map)
                // .face_patch_map(face_patch_map)
                .cosine_of_maximum_angle(cosine_of_maximum_angle)
            ;

            Mesh3 out;
            PMP::remesh_planar_patches(mesh, out, params);
            return out;
        })
        .def("interpolated_corrected_curvatures", [](
            const Mesh3& mesh,
            VertDouble& mean_curv_map,
            VertDouble& gauss_curv_map,
            VertPrincipalCurvDir& princ_curv_dir_map
        ) {
            auto params = PMP::parameters::
                vertex_mean_curvature_map(mean_curv_map)
                .vertex_Gaussian_curvature_map(mean_curv_map)
                .vertex_principal_curvatures_and_directions_map(princ_curv_dir_map)
            ;
            PMP::interpolated_corrected_curvatures(mesh, params);
        })
        .def("interpolated_corrected_curvatures", [](
            const Mesh3& mesh,
            VertDouble& mean_curv_map,
            VertDouble& gauss_curv_map,
            VertPrincipalCurvDir& princ_curv_dir_map,
            const double ball_radius
        ) {
            auto params = PMP::parameters::
                vertex_mean_curvature_map(mean_curv_map)
                .vertex_Gaussian_curvature_map(mean_curv_map)
                .vertex_principal_curvatures_and_directions_map(princ_curv_dir_map)
                .ball_radius(ball_radius)
            ;
            PMP::interpolated_corrected_curvatures(mesh, params);
        })
    ;

    py::class_<AdaptiveSizingField>(sub, "AdaptiveSizingField", py::module_local())
        .def(
            py::init([](
                const double tol,
                const std::pair<double, double>& edge_len_min_max,
                const Faces& faces,
                Mesh3& mesh
            ) {
                return AdaptiveSizingField(tol, edge_len_min_max, faces, mesh);
            }),
            py::arg("tol"),
            py::arg("edge_len_min_max"),
            py::arg("faces"),
            py::arg("mesh")
        )
        .def(
            py::init([](
                const double tol,
                const std::pair<double, double>& edge_len_min_max,
                const Faces& faces,
                Mesh3& mesh,
                double ball_radius
            ) {
                return AdaptiveSizingField(tol, edge_len_min_max, faces, mesh, PMP::parameters::ball_radius(ball_radius));
            }),
            py::arg("tol"),
            py::arg("edge_len_min_max"),
            py::arg("faces"),
            py::arg("mesh"),
            py::arg("ball_radius")
        )
    ;

    py::class_<UniformSizingField>(sub, "UniformSizingField", py::module_local())
        .def(py::init<const double, const Mesh3&>())
    ;

    // TODO expoint mesh point map
    py::class_<VertexPointMapWrapper>(sub, "VertexPointMapWrapper")
        .def(py::init<VertPoint&, VertBool&>())
    ;
}
