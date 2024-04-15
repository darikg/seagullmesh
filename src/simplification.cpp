#include "seagullmesh.hpp"
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Face_count_ratio_stop_predicate.h>
// #include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_length_stop_predicate.h>

namespace SMS = CGAL::Surface_mesh_simplification;


typedef SMS::Edge_count_stop_predicate<Mesh3>           EC;
typedef SMS::Face_count_stop_predicate<Mesh3>           FC;
typedef SMS::Edge_count_ratio_stop_predicate<Mesh3>     ECR;
typedef SMS::Face_count_stop_predicate<Mesh3>           FCR;

typedef Mesh3::Property_map<E, bool>                    EdgeBool;


template <typename T>
auto do_edge_collapse(Mesh3& mesh, const T& stop_policy, EdgeBool& edge_is_constrained) {
    auto np = CGAL::parameters::edge_is_constrained_map(edge_is_constrained);
    return SMS::edge_collapse(mesh, stop_policy, np);
}


void init_simplification(py::module &m) {
    py::module sub = m.def_submodule("simplification");
    sub
        .def("edge_collapse_edge_count", [](Mesh3& mesh, const Mesh3::size_type n, EdgeBool& edge_is_constrained) {
            return do_edge_collapse(mesh, EC(n), edge_is_constrained);
        })
        .def("edge_collapse_face_count", [](Mesh3& mesh, const Mesh3::size_type n, EdgeBool& edge_is_constrained) {
            return do_edge_collapse(mesh, FC(n), edge_is_constrained);
        })
        .def("edge_collapse_edge_count_ratio", [](Mesh3& mesh, const double ratio, EdgeBool& edge_is_constrained) {
            return do_edge_collapse(mesh, ECR(ratio), edge_is_constrained);
        })
        .def("edge_collapse_face_count_ratio", [](Mesh3& mesh, const double ratio, EdgeBool& edge_is_constrained) {
            return do_edge_collapse(mesh, FCR(ratio), edge_is_constrained);
        })
    ;
}
