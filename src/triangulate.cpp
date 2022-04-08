#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Triple<size_t, size_t, size_t> Triangle;

py::array_t<size_t> constrained_contour_pair_mesh(
    const std::vector<Point3>& p,
    const std::vector<Point3>& q,
    const std::vector<size_t>& pidx,
    const std::vector<size_t>& qidx,
    const size_t np0,
    const size_t nq0
) {
    const auto np = p.size(), nq = q.size();
    const auto n = pidx.size();
    size_t nf = 0;
    py::array_t<size_t, py::array::c_style> faces({size_t(np + nq), size_t(3)});
    auto r = faces.mutable_unchecked<2>();

    for (auto i = 0; i < (n - 1); i++) {
        const auto p0 = pidx[i], p1 = pidx[i + 1], q0 = qidx[i], q1 = qidx[i + 1];
        // p1 is always > p0 unless p1 is 0
        const auto npi = (p1 != 0) ? p1 - p0 : np - p0;
        const auto nqi = (q1 != 0) ? q1 - q0 : nq - q0;

        // Construct the border of the polygonal face to be triangulated
        // Patch has (npi + 1) P vertices, and (nqi + 1) Q vertices
        // Note that we iterate in reverse order over the q pts
        std::vector<Point3> polygon;
        polygon.reserve(npi + nqi);
        for (auto j = 0; j <= npi; j++) {
            polygon.emplace_back(p[(p0 + j) % np]);
        }
        for (auto j = nqi + 1; j-- > 0;) {
            polygon.emplace_back(q[(q0 + j) % nq]);
        }
        
        std::vector<Triangle> patch;
        patch.reserve(npi + nqi - 2);
        PMP::triangulate_hole_polyline(polygon, std::back_inserter(patch));

        // Translate the local patch back into points indices
        for (auto j = 0; j < patch.size(); j++) {
            const auto a = patch[j].first, b = patch[j].second, c = patch[j].third;
            // The q indices are a little hairy because of the reverse ordering.
            // Let v >= (npi + 1) be an index into one of the Q vertices in the patch
            //      q_patch_idx = (v - (npi + 1))  # account for the (npi + 1) P vertices
            // Because of the reverse ordering of the Q points,
            //      q_pts_idx = q0 + (nqi - q_patch_idx)
            //                = q0 + npi + nqi + 1 - v
            r(nf, 0) = (a <= npi) ? ((p0 + a) % np) + np0 : ((q0 + npi + nqi + 1 - a) % nq) + nq0;
            r(nf, 1) = (b <= npi) ? ((p0 + b) % np) + np0 : ((q0 + npi + nqi + 1 - b) % nq) + nq0;
            r(nf, 2) = (c <= npi) ? ((p0 + c) % np) + np0 : ((q0 + npi + nqi + 1 - c) % nq) + nq0;
            nf++;
        }
    }
    return faces;
}

void init_triangulate(py::module &m) {
    m.def_submodule("triangulate")
        .def("constrained_contour_pair_mesh", [](
            const py::array_t<double>& p_in, const py::array_t<double>& q_in,
            const std::vector<size_t>& pidx, const std::vector<size_t>& qidx, const size_t np0, const size_t nq0
        ) {
            std::vector<Point3> p = array_to_points_3(p_in);
            std::vector<Point3> q = array_to_points_3(q_in);
            return constrained_contour_pair_mesh(p, q, pidx, qidx, np0, nq0);
        })
    ;
}
