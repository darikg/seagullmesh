#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/poisson_surface_reconstruction.h>


typedef std::pair<Point3, Vector3> PV;

void init_poisson_reconstruct(py::module &m) {
    m.def_submodule("poisson_reconstruct")
        .def("reconstruct_surface", [](
                const py::array_t<double>& points,
                const py::array_t<double>& normals,
                const double spacing) {

            auto p = points.unchecked<2>();
            auto n = normals.unchecked<2>();
            auto np = p.shape(0);

            std::vector<PV> point_vecs;
            point_vecs.reserve(np);
            for (auto i = 0; i < np; i++) {
                point_vecs.emplace_back(std::make_pair(
                    Point3(p(i, 0), p(i, 1), p(i, 2)),
                    Vector3(n(i, 0), n(i, 1), n(i, 2))
                ));
            }

            Mesh3 mesh;

            bool success = CGAL::poisson_surface_reconstruction_delaunay(
                point_vecs.begin(),
                point_vecs.end(),
                CGAL::First_of_pair_property_map<PV>(),
                CGAL::Second_of_pair_property_map<PV>(),
                mesh,
                spacing
            );
            if (!success) {
                throw std::runtime_error("Boolean operation failed.");
            }
            return mesh;
        })
    ;
}