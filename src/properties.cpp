#include "seagullmesh.hpp"

template <typename Key, typename Val>
auto add_property_map(Mesh3& mesh, std::string name, const Val default_val) {
    typename Mesh3::Property_map<Key, Val> pmap;
    bool created;
    std::tie(pmap, created) = mesh.add_property_map<Key, Val>(name, default_val);
    if (!created) {
        throw std::runtime_error("Property map already exists");
    }
    return pmap;
}

template <typename Key, typename Val>
auto define_property_map(py::module &m, std::string name) {
    // https://stackoverflow.com/a/47749076/7519203
    using PMap = typename Mesh3::Property_map<Key, Val>;
    return py::class_<PMap>(m, name.c_str(), py::buffer_protocol(), py::dynamic_attr())
        .def("__getitem__", [](const PMap& pmap, const Key& key) {
            Val val = pmap[key];
            return val;
        })
        .def("__getitem__", [](const PMap& pmap, const std::vector<Key>& keys) {
            size_t nk = keys.size();
            py::array_t<Val, py::array::c_style> vals({int(nk)});
            auto r = vals.template mutable_unchecked<1>();

            for (size_t i = 0; i < nk; i++) {
                r(i) = pmap[keys[i]];
            }
            return vals;
        })
        .def("__setitem__", [](PMap& pmap, const Key& key, const Val val) {
            pmap[key] = val;
        })
        .def("__setitem__", [](PMap& pmap, const std::vector<Key>& keys, const Val val) {
            for (Key key : keys) {
                pmap[key] = val;
            }
        })
        .def("__setitem__", [](PMap& pmap, const std::vector<Key>& keys, const std::vector<Val>& vals) {
            size_t nk = keys.size();
            size_t nv = vals.size();
            if (nk != nv) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            for (size_t i = 0; i < nk; i++) {
                pmap[keys[i]] = vals[i];
            }
        })
    ;
}


template <typename Key, typename Val>
void define_array_3_property_map(py::module &m, std::string name) {
    using PMap = typename Mesh3::Property_map<Key, Val>;

    define_property_map<Key, Val>(m, name)
        .def("get_array", [](const PMap& pmap, const std::vector<Key>& keys) {
            const size_t nk = keys.size();
            py::array_t<double, py::array::c_style> vals({nk, size_t(3)});
            auto r = vals.mutable_unchecked<2>();

            for (auto i = 0; i < nk; i++) {
                auto val = pmap[keys[i]];
                for (auto j = 0; j < 3; j++) {
                    r(i, j) = val[j];
                }
            }
            return vals;
        })
        .def("set_array", [](PMap& pmap, const std::vector<Key>& keys, const py::array_t<double>& vals) {
            const size_t nk = keys.size();
            auto r = vals.unchecked<2>();
            if (nk != r.shape(0)) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            if (3 != r.shape(1)) {
                throw std::runtime_error("Expected an array with 3 columns");
            }
            for (auto i = 0; i < nk; i++) {
                pmap[keys[i]] = Val(r(i, 0), r(i, 1), r(i, 2));
            }
        })
    ;
}

template <typename Key, typename Val>
void define_array_2_property_map(py::module &m, std::string name) {
    using PMap = typename Mesh3::Property_map<Key, Val>;

    define_property_map<Key, Val>(m, name)
        .def("get_array", [](const PMap& pmap, const std::vector<Key>& keys) {
            const size_t nk = keys.size();
            py::array_t<double, py::array::c_style> vals({nk, size_t(2)});
            auto r = vals.mutable_unchecked<2>();

            for (auto i = 0; i < nk; i++) {
                auto val = pmap[keys[i]];
                for (auto j = 0; j < 2; j++) {
                    r(i, j) = val[j];
                }
            }
            return vals;
        })
        .def("set_array", [](PMap& pmap, const std::vector<Key>& keys, const py::array_t<double>& vals) {
            const size_t nk = keys.size();
            auto r = vals.unchecked<2>();
            if (nk != r.shape(0)) {
                throw std::runtime_error("Key and value array sizes do not match");
            }
            if (2 != r.shape(1)) {
                throw std::runtime_error("Expected an array with 2 columns");
            }
            for (auto i = 0; i < nk; i++) {
                pmap[keys[i]] = Val(r(i, 0), r(i, 1));
            }
        })
    ;
}


void init_properties(py::module &m) {
    py::module sub = m.def_submodule("properties");

    define_property_map<V, bool    >(sub, "VertBoolPropertyMap");
    define_property_map<V, int     >(sub, "VertIntPropertyMap");
    define_property_map<V, double  >(sub, "VertDoublePropertyMap");
    define_property_map<F, bool    >(sub, "FaceBoolPropertyMap");
    define_property_map<F, int     >(sub, "FaceIntPropertyMap");
    define_property_map<F, double  >(sub, "FaceDoublePropertyMap");
    define_property_map<E, bool    >(sub, "EdgeBoolPropertyMap");
    define_property_map<E, int     >(sub, "EdgeIntPropertyMap");
    define_property_map<E, double  >(sub, "EdgeDoublePropertyMap");
    define_property_map<H, bool    >(sub, "HalfedgeBoolPropertyMap");
    define_property_map<H, int     >(sub, "HalfedgeIntPropertyMap");
    define_property_map<H, double  >(sub, "HalfedgeDoublePropertyMap");

    define_array_3_property_map<V, Point3  >(sub, "VertPoint3PropertyMap");
    define_array_3_property_map<V, Vector3 >(sub, "VertVector3PropertyMap");
    define_array_2_property_map<V, Point2  >(sub, "VertPoint2PropertyMap");
    define_array_2_property_map<V, Vector2 >(sub, "VertVector2PropertyMap");

    define_array_3_property_map<F, Point3  >(sub, "FacePoint3PropertyMap");
    define_array_3_property_map<F, Vector3 >(sub, "FaceVector3PropertyMap");
    define_array_2_property_map<F, Point2  >(sub, "FacePoint2PropertyMap");
    define_array_2_property_map<F, Vector2 >(sub, "FaceVector2PropertyMap");

    define_array_3_property_map<E, Point3  >(sub, "EdgePoint3PropertyMap");
    define_array_3_property_map<E, Vector3 >(sub, "EdgeVector3PropertyMap");
    define_array_2_property_map<E, Point2  >(sub, "EdgePoint2PropertyMap");
    define_array_2_property_map<E, Vector2 >(sub, "EdgeVector2PropertyMap");

    define_array_3_property_map<H, Point3  >(sub, "HalfedgePoint3PropertyMap");
    define_array_3_property_map<H, Vector3 >(sub, "HalfedgeVector3PropertyMap");
    define_array_2_property_map<H, Point2  >(sub, "HalfedgePoint2PropertyMap");
    define_array_2_property_map<H, Vector2 >(sub, "HalfedgeVector2PropertyMap");

    sub
     .def("add_vertex_property",   &add_property_map<V, bool>)
     .def("add_vertex_property",   &add_property_map<V, int>)
     .def("add_vertex_property",   &add_property_map<V, double>)
     .def("add_vertex_property",   &add_property_map<V, Point3>)
     .def("add_vertex_property",   &add_property_map<V, Vector3>)
     .def("add_vertex_property",   &add_property_map<V, Point2>)
     .def("add_vertex_property",   &add_property_map<V, Vector2>)

     .def("add_face_property",     &add_property_map<F, bool>)
     .def("add_face_property",     &add_property_map<F, int>)
     .def("add_face_property",     &add_property_map<F, double>)
     .def("add_face_property",     &add_property_map<F, Point3>)
     .def("add_face_property",     &add_property_map<F, Vector3>)
     .def("add_face_property",     &add_property_map<F, Point2>)
     .def("add_face_property",     &add_property_map<F, Vector2>)

     .def("add_edge_property",     &add_property_map<E, bool>)
     .def("add_edge_property",     &add_property_map<E, int>)
     .def("add_edge_property",     &add_property_map<E, double>)
     .def("add_edge_property",     &add_property_map<E, Point3>)
     .def("add_edge_property",     &add_property_map<E, Vector3>)
     .def("add_edge_property",     &add_property_map<E, Point2>)
     .def("add_edge_property",     &add_property_map<E, Vector2>)

     .def("add_halfedge_property", &add_property_map<H, bool>)
     .def("add_halfedge_property", &add_property_map<H, int>)
     .def("add_halfedge_property", &add_property_map<H, double>)
     .def("add_halfedge_property", &add_property_map<H, Point3>)
     .def("add_halfedge_property", &add_property_map<H, Vector3>)
     .def("add_halfedge_property", &add_property_map<H, Point2>)
     .def("add_halfedge_property", &add_property_map<H, Vector2>)
    ;
}