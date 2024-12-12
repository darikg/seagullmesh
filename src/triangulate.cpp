#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

namespace PMP = CGAL::Polygon_mesh_processing;
typedef CGAL::Triple<size_t, size_t, size_t> Triangle;
typedef Mesh3::Property_map<V, double> VertDouble;

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


class TubeMesher {
    private:

    Mesh3& mesh;
    VertDouble& t_map;
    VertDouble& theta_map;
    FaceBool& is_cap_map;

    std::map<double, V> prev_xs;

    std::map<double, V> add_points_and_radial_edges(
            const double t,
            const py::array_t<double>& theta,
            const py::array_t<double>& pts,
    ) {
        const size_t n = pts.shape(0);
        auto r_pts = pts.unchecked<2>();
        auto r_theta = theta.unchecked<1>();
        std::map<double, V> out;

        for (size_t i = 0; i < n; ++i) {
            Point3 pt = Point3(r_pts(i, 0), r_pts(i, 1), r_pts(i, 2));
            V v = mesh.add_vertex(pt);
            t_map[v] = t;
            theta_map[v] = r_theta(i);
            out[r_theta(i)] = v;
        }
        out[2 * CGAL_PI] = out[0];  // Close the loop
        return out;
    }

    struct TriangulateCapVisitor : public PMP::PMPTriangulateFaceVisitor<Mesh3>, public PMP::PMPHolefillingVisitor<Mesh>  {
        void after_subface_created(F f) {is_cap_map[f] = true;}
    };

    public:

    TubeMesher(
            Mesh3& mesh,
            VertDouble& t_map,
            VertDouble& theta_map,
            FaceBool& is_cap_map,
        ) : mesh(mesh), t_map(t_map), theta_map(theta_map), is_cap_map(is_cap_map) {}

    void add_xs(double t, const py::array_t<double>& theta, const py::array_t<double>& pts) {
        std::map<double, V> next_xs = add_points_and_radial_edges(t, theta, pts);
        if (prev_xs.size() == 0) {
            // Must be the first xs
            prev_xs = next_xs;
            return;
        }

        std::vector<V> face;
        double theta0 = 0, theta1 = 0;
        while (theta0 < 2 * CGAL_PI) {
            // Initialize face with axial edge from prev_xs to next_xs
            face.clear();
            face.push_back(prev_xs[theta0]);
            face.push_back(next_xs[theta0]);

            // Walk forward on next_xs until we find a theta1 that is also present in prev_xs
            for (auto it1 = next_xs.upper_bound(theta0); it1 != next_xs.end(); ++it1) {
                theta1 = it1->first;
                face.push_back(it1->second);

                if (auto it0 = prev_xs.find(theta1); it0 != prev_xs.end()) {
                    // Found matching point in prev_xs  -- Walk backwards on prev_xs until we return to theta0
                    while (it0->first > theta0) {
                        face.push_back(it0->second);
                        it0--;
                    }
                    break; // Finished this face
                }
            }
            mesh.add_face(face);
            theta0 = theta1;
        }
        prev_xs = next_xs;
    }
    void close_xs(bool flip) {
        std::vector<V> face;
        face.reserve(prev_xs.size());

        for (auto const& kv : prev_xs) {
            // NB the first vertex is stored twice, once at 0 and once at 2pi
            // don't add it twice
            if (kv.first != 2 * CGAL_PI) {
                face.push_back(kv.second);
            }
        }

        if (flip){
            std::reverse(face.begin(), face.end());
        }

        F f = mesh.add_face(face);
        PMP::triangulate_face(f, mesh, PMP::parameters::visitor(TriangulateCapVisitor()));
    }
};

void init_triangulate(py::module &m) {
    py::module sub = m.def_submodule("triangulate");

    py::class_<TubeMesher>(sub, "TubeMesher")
        .def(py::init<Mesh3&, VertDouble&, VertDouble&, double, const py::array_t<double>&, const py::array_t<double>>())
        .def("add_xs", &TubeMesher::add_xs)
        .def("close_xs", &TubeMesher::close_xs)
    ;

    sub
        .def("constrained_contour_pair_mesh", [](
            const py::array_t<double>& p_in, const py::array_t<double>& q_in,
            const std::vector<size_t>& pidx, const std::vector<size_t>& qidx, const size_t np0, const size_t nq0
        ) {
            std::vector<Point3> p = array_to_points_3(p_in);
            std::vector<Point3> q = array_to_points_3(q_in);
            return constrained_contour_pair_mesh(p, q, pidx, qidx, np0, nq0);
        })
        .def("triangulate_faces", [](Mesh3& mesh, const std::vector<F>& faces) {
            PMP::triangulate_faces(faces, mesh);
        })
        .def("reverse_face_orientations", [](Mesh3& mesh, const std::vector<F>& faces) {
            PMP::reverse_face_orientations(faces, mesh);
        })
        .def("does_bound_a_volume", [](const Mesh3& mesh) {
            return PMP::does_bound_a_volume(mesh);
        })
        .def("is_outward_oriented", [](const Mesh3& mesh) {
            return PMP::is_outward_oriented(mesh);
        })
    ;
}
