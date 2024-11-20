#include "seagullmesh.hpp"
#include "util.hpp"

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/IO/OFF.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Bbox_3 BBox3;

template<typename T>
void define_simple_type_3(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str(), py::module_local())
        .def(py::init<double, double, double>())
        .def("__eq__", [](const T& self, const T& other) {return self == other;})
    ;
}

template<typename T>
void define_simple_type_2(py::module &m, std::string name) {
    py::class_<T>(m, name.c_str(), py::module_local())
        .def(py::init<double, double>())
        .def("__eq__", [](const T& self, const T& other) {return self == other;})
    ;
}

// Used for defining Vertex/Vertices, Face/Faces, etc.
template<typename Idx>
void define_indices(py::module &m, std::string idx_name, std::string idxs_name) {
    using size_type = typename Idx::size_type;
    using Idxs = typename std::vector<Idx>;

    py::class_<Idx>(m, idx_name.c_str())
        .def(py::init<size_type>())
        .def("to_int", [](const Idx& idx) {return size_type(idx);})
    ;
    py::class_<Idxs>(m, idxs_name.c_str())
        .def("__len__", [](const Idxs& idxs) { return idxs.size(); })
        .def("__iter__", [](Idxs& idxs) {
                return py::make_iterator(idxs.begin(), idxs.end());
            }, py::keep_alive<0, 1>()  /* Keep vector alive while iterator is used */
        )
        .def("__getitem__", [](const Idxs& idxs, size_t i) {
            if (i >= idxs.size()) {
                throw py::index_error();
            }
            return idxs[i];
        })
        .def("__getitem__", [](const Idxs& idxs, py::slice slice) {
            py::ssize_t start, stop, step, slicelength;
            if (!slice.compute(idxs.size(), &start, &stop, &step, &slicelength)) {
                throw py::error_already_set();
            }
            Idxs out;
            out.reserve(slicelength);
            for (int i = 0; i < slicelength; ++i) {
                out.emplace_back(idxs[start]);
                start += step;
            }
            return out;
        })
        .def("__getitem__", [](const Idxs& idxs, const py::array_t<size_type>& sub) {
            if (sub.ndim() != 1) {
                throw py::index_error();
            }
            py::ssize_t n = sub.size();
            Idxs out;
            out.reserve(n);

            auto r = sub.template unchecked<1>();
            for (int i = 0; i < n; ++i) {
                out.emplace_back(idxs.at(r(i)));
            }
            return out;
        })
        .def("__getitem__", [](const Idxs& idxs, const py::array_t<bool>& sub) {
            py::ssize_t n = sub.size();
            if (sub.ndim() != 1 || n != idxs.size()) {
                throw py::index_error();
            }
            Idxs out;
            auto r = sub.template unchecked<1>();
            for (int i = 0; i < n; ++i) {
                if ( r(i) ) {
                    out.emplace_back(idxs[i]);
                }
            }
            return out;
        })
        .def("__eq__", [](const Idxs& idxs, const Idxs& other) {
            if (idxs.size() != other.size()) {
                return false;
            }
            for (int i = 0; i < idxs.size(); ++i) {
                if ( idxs[i] != other[i] ) {
                    return false;
                }
            }
            return true;
        })
        .def("to_ints", [](const Idxs& idxs) {
            const py::ssize_t n = idxs.size();
            py::template array_t<size_type> out({n});
            auto r = out.template mutable_unchecked<1>();
            for (int i = 0; i < n; ++i) {
                r(i) = size_type(idxs[i]);
            }
            return out;
        })
        .def("from_ints", [](const py::array_t<size_type>& idxs) {
            if (idxs.ndim() != 1) {
                throw py::index_error();
            }
            py::ssize_t n = idxs.size();
            Idxs out;
            out.reserve(n);
            auto r = idxs.template unchecked<1>();
            for (int i = 0; i < n; ++i) {
                out.emplace_back(V(r(i)));
            }
            return out;
        })
    ;
}

template<typename Idx, typename IdxRange>
std::vector<Idx> indices_from_range(typename Idx::size_type n, const IdxRange idx_range) {
    std::vector<Idx> idxs;
    idxs.reserve(n);
    for (Idx idx : idx_range) {
        idxs.emplace_back(idx);
    }
    return idxs;
}


template<size_t N, typename Idx, typename Val>
py::array_t<double> map_indices_to_array(const std::vector<Idx>& idxs, const std::function<Val (Idx)> fn) {
    const size_t n_idxs = idxs.size();
    py::array_t<double> vals({n_idxs, N});
    auto r = vals.mutable_unchecked<2>();
    for (auto i = 0; i < n_idxs; ++i) {
        Val val = fn(idxs[i]);
        for (auto j = 0; j < N; ++j) {
            r(i, j) = CGAL::to_double(val[j]);
        }
    }
    return vals;
}


void init_mesh(py::module &m) {
    py::module sub = m.def_submodule("mesh");

    define_simple_type_2<Point2>(sub, "Point2");
    define_simple_type_3<Point3>(sub, "Point3");
    define_simple_type_2<Vector2>(sub, "Vector2");
    define_simple_type_3<Vector3>(sub, "Vector3");

    define_indices<V>(sub, "Vertex", "Vertices");
    define_indices<F>(sub, "Face", "Faces");
    define_indices<E>(sub, "Edge", "Edges");
    define_indices<H>(sub, "Halfedge", "Halfedges");

    sub.def("polygon_soup_to_mesh3", [](
                py::array_t<double> &points,
                std::vector<std::vector<size_t>>& faces,
                const bool orient
        ) {
            Mesh3 mesh;
            std::vector<Point3> vertices = array_to_points_3(points);

            if (orient) {
                bool success = PMP::orient_polygon_soup(vertices, faces);
                if (!success) {
                    throw std::runtime_error("Polygon orientation failed");
                }
            }
            PMP::polygon_soup_to_polygon_mesh(vertices, faces, mesh);
            return mesh;
        })
        .def("load_mesh_from_file", [](const std::string filename) {
            Mesh3 mesh;
            if(!CGAL::IO::read_polygon_mesh(filename, mesh)) {
                throw std::runtime_error("Failed to load mesh");
            }
            return mesh;
        })
    ;

    py::class_<BBox3>(sub, "BoundingBox3", py::module_local())
        .def(py::init<double, double, double, double, double, double>())
        .def_property_readonly("x_min", &BBox3::xmin)
        .def_property_readonly("x_max", &BBox3::xmax)
        .def_property_readonly("y_min", &BBox3::ymin)
        .def_property_readonly("y_max", &BBox3::ymax)
        .def_property_readonly("x_min", &BBox3::zmin)
        .def_property_readonly("z_max", &BBox3::zmax)
        .def("diagonal", [](const BBox3& bbox) {
            return std::sqrt(
                CGAL::square(bbox.xmax() - bbox.xmin()) +
                CGAL::square(bbox.ymax() - bbox.ymin()) +
                CGAL::square(bbox.zmax() - bbox.zmin())
            );
        })
    ;

    py::class_<Mesh3>(sub, "Mesh3")
        .def(py::init<>())

        .def_property_readonly("has_garbage", [](const Mesh3& mesh) {return mesh.has_garbage();})
        .def("collect_garbage", [](Mesh3& mesh) {mesh.collect_garbage();})
        
        .def_property_readonly("is_valid", [](const Mesh3& mesh) { return mesh.is_valid(false); })
        .def_property_readonly("n_vertices", [](const Mesh3& mesh) { return mesh.number_of_vertices(); })
        .def_property_readonly("n_faces", [](const Mesh3& mesh) { return mesh.number_of_faces(); })
        .def_property_readonly("n_edges", [](const Mesh3& mesh) { return mesh.number_of_edges(); })
        .def_property_readonly("n_halfedges", [](const Mesh3& mesh) { return mesh.number_of_halfedges(); })
        .def_property_readonly("points", [](const Mesh3& mesh) { return mesh.points(); })

        .def_property_readonly("vertices", [](const Mesh3& mesh) {
            return indices_from_range<V, Mesh3::Vertex_range>(mesh.number_of_vertices(), mesh.vertices());
        })
        .def_property_readonly("faces", [](const Mesh3& mesh) {
            return indices_from_range<F, Mesh3::Face_range>(mesh.number_of_faces(), mesh.faces());
        })
        .def_property_readonly("edges", [](const Mesh3& mesh) {
            return indices_from_range<E, Mesh3::Edge_range>(mesh.number_of_edges(), mesh.edges());
        })
        .def_property_readonly("halfedges", [](const Mesh3& mesh) {
            return indices_from_range<H, Mesh3::Halfedge_range>(mesh.number_of_halfedges(), mesh.halfedges());
        })

        .def("bounding_box", [](const Mesh3& mesh) {
            return PMP::bbox(mesh);
        })

        .def("edge_vertices", [](const Mesh3& mesh, const std::vector<E>& edges) {
            std::map<V, size_t> vert_idxs;
            size_t vi = 0;
            for (V v : mesh.vertices()) {
                vert_idxs[v] = vi;
                vi++;
            }

            const size_t ne = edges.size();
            py::array_t<size_t, py::array::c_style> verts({ne, size_t(2)});
            auto r = verts.mutable_unchecked<2>();
            for (auto i = 0; i < ne; i++) {
                for (auto j = 0; j < 2; j++) {
                    r(i, j) = vert_idxs[mesh.vertex(edges[i], j)];
                }
            }

            return verts;
        })
        .def("vertices_to_faces", [](const Mesh3& mesh, const std::vector<V>& verts) {
            std::set<F> faces;
            for (V v : verts) {
                for (F f : faces_around_target(mesh.halfedge(v), mesh)) {
                    if (f != mesh.null_face()) {
                        faces.insert(f);
                    }
                }
            }
            return std::vector<F>(faces.begin(), faces.end());
        })
        .def("vertices_to_edges", [](const Mesh3& mesh, const std::vector<V>& verts) {
            std::set<E> edges;
            for (V v : verts) {
                for (H h : halfedges_around_source(v, mesh)) {
                    edges.insert(mesh.edge(h));
                }
            }
            return std::vector<E>(edges.begin(), edges.end());
        })
        .def("vertex_degrees", [](const Mesh3& mesh, const std::vector<V>& verts) {
            size_t n = verts.size();
            py::array_t<Mesh3::size_type> degrees({py::ssize_t(n)});
            auto r = degrees.mutable_unchecked<1>();
            for (size_t i = 0; i < n; ++i) {
                r(i) = mesh.degree(verts[i]);
            }
            return degrees;
        })
        .def("to_polygon_soup", [](const Mesh3& mesh) {
            std::vector<Point3> verts;
            std::vector<std::vector<size_t>> faces;
            PMP::polygon_mesh_to_polygon_soup(mesh, verts, faces);
            auto points = points_to_array(verts);

            // Convert vector<vector<size_t>> to array
            const size_t nf = mesh.number_of_faces();
            py::array_t<size_t, py::array::c_style> faces_out({nf, size_t(3)});
            auto rf = faces_out.mutable_unchecked<2>();
            for (size_t i = 0; i < nf; i++) {
                for (size_t j = 0; j < 3; j++) {
                    rf(i, j) = faces[i][j];
                }
            }
            return std::make_tuple(points, faces_out);
        })
        .def("face_normals", [](const Mesh3& mesh, const std::vector<F>& faces) {
            return map_indices_to_array<3, F, Vector3>(
                faces, [&mesh](F f) {return PMP::compute_face_normal(f, mesh);}
            );
        })
        .def("vertex_normals", [](const Mesh3& mesh, const std::vector<V>& verts) {
            return map_indices_to_array<3, V, Vector3>(
                verts, [&mesh](V v) {return PMP::compute_vertex_normal(v, mesh);}
            );
        })
        .def("volume", [](const Mesh3& mesh) {return PMP::volume(mesh);})
        .def("area", [](const Mesh3& mesh) {return PMP::area(mesh);})
        .def("estimate_geodesic_distances", [](const Mesh3& mesh, Mesh3::Property_map<V, double>& distances, V source) {
            CGAL::Heat_method_3::estimate_geodesic_distances(mesh, distances, source);
        })
        .def("estimate_geodesic_distances", [](
                const Mesh3& mesh, Mesh3::Property_map<V, double>& distances, const std::vector<V>& sources) {
            CGAL::Heat_method_3::estimate_geodesic_distances(mesh, distances, sources);
        })

        .def("write_ply", [](Mesh3& mesh, std::string file) {
            std::ofstream out(file, std::ios::binary);
            CGAL::IO::set_binary_mode(out);
            bool success = CGAL::IO::write_PLY(out, mesh, "");
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })
        .def("write_off", [](Mesh3& mesh, std::string file) {
            std::ofstream out(file);
            bool success = CGAL::IO::write_OFF(out, mesh);
            if (!success) {
                throw std::runtime_error("writing failed");
            }
        })
    ;
}
