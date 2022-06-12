#include "seagullmesh.hpp"

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/fair.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/smooth_shape.h>
#include <CGAL/Polygon_mesh_processing/refine.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef std::vector<V>                     Verts;
typedef std::vector<F>                     Faces;
typedef Mesh3::Property_map<V, Point3>     VertPoint;
typedef Mesh3::Property_map<V, bool>       VertBool;


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


void init_meshing(py::module &m) {
    m.def_submodule("meshing")
        .def("remesh", [](Mesh3& mesh, const Faces& faces, double target_edge_length, unsigned int n_iter, bool protect_constraints) {
            auto params = PMP::parameters::number_of_iterations(n_iter).protect_constraints(protect_constraints);
            PMP::isotropic_remeshing(faces, target_edge_length, mesh, params);
        })
        .def("remesh", [](Mesh3& mesh, const Faces& faces, double target_edge_length, unsigned int n_iter,
                        const bool protect_constraints, VertBool& touched) {

            auto points = mesh.points();
            VertexPointMapWrapper point_map = VertexPointMapWrapper(points, touched);
            auto params = PMP::parameters::number_of_iterations(n_iter)
                .vertex_point_map(point_map)
                .protect_constraints(protect_constraints)
            ;

            PMP::isotropic_remeshing(faces, target_edge_length, mesh, params);
        })
        .def("fair", [](Mesh3& mesh, const Verts& verts, const unsigned int fairing_continuity) {
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
        .def("smooth_mesh", [](Mesh3& mesh, const std::vector<F>& faces, unsigned int n_iter, bool use_safety_constraints) {
            auto params = PMP::parameters::number_of_iterations(n_iter).use_safety_constraints(use_safety_constraints);
            PMP::smooth_mesh(faces, mesh, params);
        })
        .def("smooth_shape", [](Mesh3& mesh, const std::vector<F>& faces, const double time, unsigned int n_iter) {
            auto params = PMP::parameters::number_of_iterations(n_iter);
            PMP::smooth_shape(faces, mesh, time, params);
        })
        .def("does_self_intersect", [](Mesh3& mesh) {
            return PMP::does_self_intersect(mesh);
        })
    ;
}
