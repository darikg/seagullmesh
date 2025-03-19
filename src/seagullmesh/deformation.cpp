#include "seagullmesh.hpp"
#include "util.hpp"

// typedef CGAL::Surface_mesh_deformation<Mesh3> SurfaceMeshDeformation;
//typedef CGAL::Surface_mesh_deformation<Mesh3, CGAL::Default, CGAL::Default, CGAL::SRE_ARAP>         SMD_SreArap;
//typedef CGAL::Surface_mesh_deformation<Mesh3, CGAL::Default, CGAL::Default, CGAL::ORIGINAL_ARAP>    SMD_OriginalArap;
//typedef CGAL::Surface_mesh_deformation<Mesh3, CGAL::Default, CGAL::Default, CGAL::SPOKES_AND_RIMS>  SMD_SpokesAndRims;


template <CGAL::Deformation_algorithm_tag Tag>
auto define_deformation(py::module &m, Tag tag, std::string name) {
    using TDeform = typename CGAL::Surface_mesh_deformation<Mesh3, CGAL::Default, CGAL::Default, tag>;
    return py::class_<TDeform>(sub, name.c_str())
        .def(py::init<Mesh&>())
        .def("add_control_vertices", [](TDeform& deformation, const Vertices& vertices) {
            vertices.apply([&] (V v) { deformation.insert_control_vertex(v); });
        })
        .def("add_roi_vertices", [](TDeform& deformation, const Vertices& vertices) {
            vertices.apply([&] (V v) { deformation.insert_roi_vertex(v); });
        })
        .def("preprocess", &TDeform::preprocess)  // returns success
        .def("deform", [](TDeform& deformation, unsigned int iterations, double tolerance) {
            deformation.deform(iterations, tolerance);
        })
        .def("set_target_positions", [](
            TDeform& deformation,
            const Vertices& vertices,
            const py::array_t<double>& positions
        ) {
            auto r = positions.unchecked<2>();
            vertices.apply([&] (size_t i, V v) {
                deformation.set_target_position(v, Point3(r(i, 0), r(i, 1), r(i, 2)));
            });
        })
    ;
}


void init_connected(py::module &m) {
    py::module sub = m.def_submodule("deformation");
    define_deformation(sub, CGAL::SRE_ARAP, "SurfaceMeshDeformation_SreArap");
    define_deformation(sub, CGAL::ORIGINAL_ARAP, "SurfaceMeshDeformation_OriginalArap");
    define_deformation(sub, CGAL::SPOKES_AND_RIMS, "SurfaceMeshDeformation_SpokesAndRims");
}

