#include "seagullmesh.hpp"
#include "util.hpp"

#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Surface_mesh.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>

typedef CGAL::Mean_curvature_flow_skeletonization<Mesh3>    Skeletonization;
typedef Skeletonization::Skeleton                           Skeleton;
typedef Skeleton::vertex_descriptor                         SkelVertex;
typedef Skeleton::edge_descriptor                           SkelEdge;



void init_skeletonization(py::module &m) {
    py::module sub = m.def_submodule("skeletonization");

    /*
    The graph type representing the skeleton. The vertex property Vmap is a struct with a member point of type
    Traits::Point_3 and a member vertices of type std::vector<boost::graph_traits<TriangleMesh>::vertex_descriptor>.
    */
    py::class_<Skeleton>(sub, "Skeleton", py::module_local())
        .def_property_readonly("edges", [](const Skeleton& skeleton) {
            const auto ne = boost::num_edges(skeleton);
            py::array_t<size_t, py::array::c_style> edge_array({ne, size_t(2)});
            auto r = edge_array.mutable_unchecked<2>();
            size_t ei = 0;
            for (SkelEdge e : CGAL::make_range(edges(skeleton))) {
                r(ei, 0) = source(e, skeleton);
                r(ei, 1) = target(e, skeleton);
                ei++;
            }
            return edge_array;
        })
        .def_property_readonly("points", [](const Skeleton& skeleton) {
            const auto nv = boost::num_vertices(skeleton);
            py::array_t<double, py::array::c_style> points({nv, size_t(3)});
            auto r = points.mutable_unchecked<2>();
            size_t vi = 0;
            for (SkelVertex v : CGAL::make_range(vertices(skeleton))) {
                const Point3& p = skeleton[v].point;
                for (auto j = 0; j < 3; j++) {
                    r(vi, j) = CGAL::to_double(p[j]);
                }
                vi++;
            }
            return points;
        })
        .def_property_readonly("vertex_map", [](const Skeleton& skeleton) {
            std::map<SkelVertex, std::vector<V>> vert_map;
            for (SkelVertex v_skel : CGAL::make_range(vertices(skeleton))) {
                vert_map[v_skel] = skeleton[v_skel].vertices;
            }
            return vert_map;
        })
    ;

    sub.def("extract_mean_curvature_flow_skeleton", [](const Mesh3& mesh) {
        Skeleton skeleton;
        CGAL::extract_mean_curvature_flow_skeleton(mesh, skeleton);
        return skeleton;
    });
}