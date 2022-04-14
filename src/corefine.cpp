#include "seagullmesh.hpp"
#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef Mesh3::Property_map<V, int>        VertexIndex;
typedef Mesh3::Property_map<F, int>        FaceIndex;
typedef Mesh3::Property_map<E, bool>       EdgeConstrainedMap;


struct CorefinementVertexTracker : public PMP::Corefinement::Default_visitor<Mesh3> {
    // Used for tracking for refinement indices
    // CGAL's corefine only uses a visitor for the first mesh, so we need the references to both
    // here to tell which is which
    Mesh3& mesh1;
    Mesh3& mesh2;
    VertexIndex& vert_ids1;
    VertexIndex& vert_ids2;

    CorefinementVertexTracker(
        Mesh3& m1, Mesh3& m2, VertexIndex& v1, VertexIndex& v2
    ) : mesh1(m1), mesh2(m2), vert_ids1(v1), vert_ids2(v2) {}

    void new_vertex_added(size_t i_id, V v, const Mesh3& mesh) {
        // Called when a new vertex is added in the mesh
        // (either an edge split or a vertex inserted in the interior of a face).
        // i_id is the intersection point id reported in new_node_added.
        // For each mesh, a vertex with a given id will be reported exactly once,
        // except if it is already an existing vertex.
        if (&mesh == &mesh1) {
            vert_ids1[v] = int(i_id);
        } else {
            vert_ids2[v] = int(i_id);
        }
    }
};

struct CorefinementVertexFaceTracker : public PMP::Corefinement::Default_visitor<Mesh3> {
//https://github.com/CGAL/cgal/blob/master/Polygon_mesh_processing/examples/Polygon_mesh_processing/corefinement_mesh_union_with_attributes.cpp

    // Used for tracking for refinement indices
    // CGAL's corefine only uses a visitor for the first mesh, so we need the references to both
    // here to tell which is which
    Mesh3& mesh1;
    Mesh3& mesh2;
    VertexIndex& vert_ids1;
    VertexIndex& vert_ids2;
    FaceIndex& face_ids1;
    FaceIndex& face_ids2;
    int face_id;


    CorefinementVertexFaceTracker(
        Mesh3& m1, Mesh3& m2, VertexIndex& v1, VertexIndex& v2, FaceIndex& f1, FaceIndex&f2
    ) : mesh1(m1), mesh2(m2), vert_ids1(v1), vert_ids2(v2), face_ids1(f1), face_ids2(f2), face_id(-1) {}

    void new_vertex_added(size_t i_id, V v, const Mesh3& mesh) {
        if (&mesh == &mesh1) {
            vert_ids1[v] = int(i_id);
        } else {
            vert_ids2[v] = int(i_id);
        }
    }
    void before_subface_creations(F f_split, Mesh3& mesh) {
        if (&mesh == &mesh1) {
            face_id = face_ids1[f_split];
        } else {
            face_id = face_ids2[f_split];
        }
    }
    void after_subface_created(F f_new, Mesh3& mesh) {
        if (&mesh == &mesh1) {
            face_ids1[f_new] = face_id;
        } else {
            face_ids2[f_new] = face_id;
        }
    }
    void after_face_copy(F f_src, Mesh3& mesh_src, F f_tgt, Mesh3& mesh_tgt) {
        face_ids1[f_tgt] = face_ids2[f_src];
    }
};


void init_corefine(py::module &m) {
    py::module sub = m.def_submodule("corefine");

    py::class_<CorefinementVertexTracker>(sub, "CorefinementVertexTracker")
        .def(py::init<Mesh3&, Mesh3&, VertexIndex&, VertexIndex&>())
    ;

    py::class_<CorefinementVertexFaceTracker>(sub, "CorefinementVertexFaceTracker")
        .def(py::init<Mesh3&, Mesh3&, VertexIndex&, VertexIndex&, FaceIndex&, FaceIndex&>())
    ;


    sub.def("corefine", [](Mesh3& mesh1, Mesh3& mesh2){
        PMP::corefine(mesh1, mesh2);
    })
    .def("corefine", [](
            Mesh3& mesh1, Mesh3& mesh2, 
            EdgeConstrainedMap& ecm1, EdgeConstrainedMap& ecm2,
            CorefinementVertexTracker& tracker) {

        auto params1 = PMP::parameters::visitor(tracker).edge_is_constrained_map(ecm1);
        auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
        PMP::corefine(mesh1, mesh2, params1, params2);
    })
    .def("corefine", [](
            Mesh3& mesh1, Mesh3& mesh2,
            EdgeConstrainedMap& ecm1, EdgeConstrainedMap& ecm2,
            CorefinementVertexFaceTracker& tracker) {

        auto params1 = PMP::parameters::visitor(tracker).edge_is_constrained_map(ecm1);
        auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
        PMP::corefine(mesh1, mesh2, params1, params2);
    })
    .def("difference", [](Mesh3& mesh1, Mesh3& mesh2, Mesh3& out) {
        bool success = PMP::corefine_and_compute_difference(mesh1, mesh2, out);
        if (!success) {
            throw std::runtime_error("Boolean operation failed.");
        }
    })
    .def("union", [](Mesh3& mesh1, Mesh3& mesh2, Mesh3& out) {
        bool success = PMP::corefine_and_compute_union(mesh1, mesh2, out);
        if (!success) {
            throw std::runtime_error("Boolean operation failed.");
        }
    })
    .def("intersection", [](Mesh3& mesh1, Mesh3& mesh2, Mesh3& out) {
        bool success = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(mesh1, mesh2, out);
        if (!success) {
            throw std::runtime_error("Boolean operation failed.");
        }
    })
    .def("union", [](
            Mesh3& mesh1, Mesh3& mesh2,
            EdgeConstrainedMap& ecm1, EdgeConstrainedMap& ecm2,
            CorefinementVertexTracker& tracker) {

        auto params1 = PMP::parameters::visitor(tracker).edge_is_constrained_map(ecm1);
        auto params2 = PMP::parameters::edge_is_constrained_map(ecm2);
        bool success = PMP::corefine_and_compute_union(mesh1, mesh2, mesh1, params1, params2);
        if (!success) {
            throw std::runtime_error("Boolean operation failed.");
        }
    })
    ;
}