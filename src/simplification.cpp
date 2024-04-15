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
//    py::class_<EC>(sub, "EdgeCount").def(py::init<const E>());
//    py::class_<FC>(sub, "FaceCount").def(py::init<const F>());
//    py::class_<ECR>(sub, "EdgeCountRatio").def(py::init<const double>());
//    py::class_<FCR>(sub, "FaceCountRatio").def(py::init<const double>());

    sub
        .def("edge_count", [](const Mesh3::size_type n) {return EC(n);})
        .def("face_count", [](const Mesh3::size_type n) {return FC(n);})
        .def("edge_count_ratio", [](const double ratio) {return ECR(ratio);})
        .def("face_count_ratio", [](const double ratio) {return FCR(ratio);})
//        .def("edge_collapse", [](Mesh3& mesh, const EC& stop_policy, EdgeBool& edge_is_constrained) {
//            return do_edge_collapse(mesh, stop_policy, edge_is_constrained);
//        })
//        .def("edge_collapse", [](Mesh3& mesh, const FC& stop_policy, EdgeBool& edge_is_constrained) {
//            return do_edge_collapse(mesh, stop_policy, edge_is_constrained);
//        })
//        .def("edge_collapse", [](Mesh3& mesh, const ECR& stop_policy, EdgeBool& edge_is_constrained) {
//            return do_edge_collapse(mesh, stop_policy, edge_is_constrained);
//        })
//        .def("edge_collapse", [](Mesh3& mesh, const FCR& stop_policy, EdgeBool& edge_is_constrained) {
//            return do_edge_collapse(mesh, stop_policy, edge_is_constrained);
//        })
    ;
}
