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


void init_simplification(py::module &m) {
    py::module sub = m.def_submodule("simplification");
    py::class_<EC>(sub, "EdgeCount").def(py::init<const E>());
    py::class_<FC>(sub, "FaceCount").def(py::init<const F>());
    py::class_<ECR>(sub, "EdgeCountRatio").def(py::init<const double>());
    py::class_<FCR>(sub, "FaceCountRatio").def(py::init<const double>());

    sub
        .def("edge_collapse", [](Mesh3& mesh, EC stop) {SMS::edge_collapse(mesh, stop);})
    ;
}
