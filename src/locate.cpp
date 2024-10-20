#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/Polygon_mesh_processing/locate.h>
#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_tree.h>

// TODO AABB_traits is deprecated, use AABB_traits_3
#include <CGAL/AABB_traits.h>
//#include <CGAL/AABB_traits_3.h>

typedef std::vector<F>         Faces3;
typedef std::vector<V>         Verts3;

namespace PMP = CGAL::Polygon_mesh_processing;
typedef PMP::Barycentric_coordinates<Kernel::FT>    Barycentric_coordinates;
typedef PMP::Face_location<Mesh3, Kernel::FT>       FaceLocation;

typedef Mesh3::Property_map<V, Point3>             VertPoints3;
typedef Mesh3::Property_map<V, Point2>             VertPoints2;

typedef typename CGAL::AABB_face_graph_triangle_primitive<Mesh3>    AABB_primitive3;
typedef typename CGAL::AABB_traits<Kernel, AABB_primitive3>         AABB_traits3;
typedef typename CGAL::AABB_tree<AABB_traits3>                      AABB_Tree3;

template<typename Point, typename VPM>
auto construct_points(
        const Mesh3& mesh,
        const std::vector<F>& faces,
        const py::array_t<double>& bary_coords,
        const VPM& vertex_point_map
) {
    size_t nf = faces.size();
    auto rbc = bary_coords.unchecked<2>();
    size_t nb = rbc.shape(0);
    if (nf != nb) {
        throw std::runtime_error("number of faces doesn't match number of points");
    }

    std::vector<Point> points;
    points.reserve(nf);
    auto params = CGAL::parameters::vertex_point_map(vertex_point_map);

    for (auto i = 0; i < nf; i++) {
        Barycentric_coordinates bc = {rbc(i, 0), rbc(i, 1), rbc(i, 2)};
        FaceLocation loc = {faces[i], bc};
        auto pt = PMP::construct_point(loc, mesh, params);
        points.emplace_back(pt);
    }
    return points_to_array(points);
}

struct Point2_to_Point3 {
    // https://stackoverflow.com/questions/66308313/2d-aabbtree-in-cgal-with-custom-property-map
    // https://stackoverflow.com/questions/24317345/cgal-using-locate-to-find-cell-on-triangulation-surface
    using key_type = V;
    using value_type = Point3;
    using reference = Point3;
    using category = boost::readable_property_map_tag;

    const VertPoints2* points;

    Point2_to_Point3() : points(nullptr) {}
    Point2_to_Point3(const VertPoints2& points) : points(&points) {}

    friend Point3 get(const Point2_to_Point3 &map, V v) {
        auto p = map.points->operator[](v);
        return {p[0], p[1], 0};
    }
};

typedef typename CGAL::AABB_face_graph_triangle_primitive<Mesh3, Point2_to_Point3>      AABB_primitive2;
typedef typename CGAL::AABB_traits<Kernel, AABB_primitive2>                             AABB_traits2;
typedef typename CGAL::AABB_tree<AABB_traits2>                                          AABB_Tree2;

template<typename AABB_Tree, typename Point, typename VPM>
auto locate_points(
        const Mesh3& mesh,
        const AABB_Tree& tree,
        const std::vector<Point>& points,
        const VPM& vertex_point_map
) {
    size_t np = points.size();
    std::vector<F> faces;
    faces.reserve(np);
    py::array_t<double, py::array::c_style> bary_coords({np, size_t(3)});
    auto params = CGAL::parameters::vertex_point_map(vertex_point_map);

    auto rbc = bary_coords.mutable_unchecked<2>();

    for (auto i = 0; i < np; i++) {
        FaceLocation loc = PMP::locate_with_AABB_tree(points[i], tree, mesh, params);
        faces.emplace_back(loc.first);

        for (size_t j = 0; j < 3; j++) {
            rbc(i, j) = loc.second[j];
        }
    }

    return std::make_tuple(faces, bary_coords);
}

void init_locate(py::module &m) {
    py::module sub = m.def_submodule("locate");

    py::class_<AABB_Tree3>(sub, "AABB_Tree3");

    sub.def("aabb_tree", [](const Mesh3& mesh) {
            AABB_Tree3 tree;
            PMP::build_AABB_tree(mesh, tree);
            return tree;
        })
        .def("aabb_tree", [](const Mesh3& mesh, const VertPoints3& point_map) {
            AABB_Tree3 tree;
            auto params = CGAL::parameters::vertex_point_map(point_map);
            PMP::build_AABB_tree(mesh, tree, params);
            return tree;
        })
        .def("aabb_tree", [](const Mesh3& mesh, const VertPoints2& points2) {
            AABB_Tree2 tree;
            Point2_to_Point3 points3(points2);

            auto params = CGAL::parameters::vertex_point_map(points3);
            PMP::build_AABB_tree(mesh, tree, params);
            return tree;
        })
        .def("locate_points", [](const Mesh3& mesh, const AABB_Tree3& tree, const py::array_t<double>& points) {
            auto pts = array_to_points_3(points);
            return locate_points(mesh, tree, pts, mesh.points());
        })
        .def("locate_points", [](
                const Mesh3& mesh,
                const AABB_Tree3& tree,
                const py::array_t<double>& points,
                const VertPoints3& vertex_point_map
            ) {
            auto pts = array_to_points_3(points);
            return locate_points(mesh, tree, pts, vertex_point_map);
        })
        .def("locate_points", [](
                const Mesh3& mesh,
                const AABB_Tree2& tree,
                const py::array_t<double>& points,
                const VertPoints2& vertex_point_map
        ) {
            auto pts = array_to_points_2(points);
            return locate_points(mesh, tree, pts, vertex_point_map);
        })
        .def("construct_points", [](
                const Mesh3& mesh,
                const std::vector<F>& faces,
                const py::array_t<double>& bary_coords,
                const VertPoints3& vertex_point_map
        ){
            return construct_points<Point3, VertPoints3>(mesh, faces, bary_coords, mesh.points());
        })
        .def("construct_points", [](
                const Mesh3& mesh,
                const std::vector<F>& faces,
                const py::array_t<double>& bary_coords,
                const VertPoints3& vertex_point_map
        ){
            return construct_points<Point3, VertPoints3>(mesh, faces, bary_coords, vertex_point_map);
        })
        .def("construct_points", [](
                const Mesh3& mesh,
                const std::vector<F>& faces,
                const py::array_t<double>& bary_coords,
                const VertPoints2& vertex_point_map
        ){
            return construct_points<Point2, VertPoints2>(mesh, faces, bary_coords, vertex_point_map);
        })
        .def("shortest_path", [](
                const Mesh3& mesh,
                const F src_face, const std::vector<double>& src_bc,
                const F tgt_face, const std::vector<double>& tgt_bc) {

            using SPTraits = CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh3>;
            using ShortestPath = CGAL::Surface_mesh_shortest_path<SPTraits>;

            Barycentric_coordinates src_bc_ = {src_bc[0], src_bc[1], src_bc[2]};
            Barycentric_coordinates tgt_bc_ = {tgt_bc[0], tgt_bc[1], tgt_bc[2]};

            ShortestPath shortest_path(mesh);
            shortest_path.add_source_point(src_face, src_bc_);
            std::vector<Point3> points;
            shortest_path.shortest_path_points_to_source_points(tgt_face, tgt_bc_, std::back_inserter(points));

            return points_to_array(points);
        })
    ;
}
