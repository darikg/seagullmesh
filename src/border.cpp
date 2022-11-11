#include "seagullmesh.hpp"
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/iterator/function_output_iterator.hpp>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef Mesh3::Property_map<V, bool>       VertBool;
typedef Mesh3::Property_map<E, bool>       EdgeBool;


struct touch_border_vertices {
  const Mesh3& mesh;
  VertBool& verts;

  touch_border_vertices (const Mesh3& m, VertBool& v) : mesh(m), verts(v) {}

  void operator()(const H& h) const {
    verts[mesh.source(h)] = true;
    verts[mesh.target(h)] = true;
  }

};


void init_border(py::module &m) {
    m.def_submodule("border")
        .def("label_border_vertices", [](const Mesh3& mesh, VertBool& vert_is_border) {
            auto output_iter = boost::make_function_output_iterator(touch_border_vertices(mesh, vert_is_border));
            PMP::border_halfedges(faces(mesh), mesh, output_iter);
        })
    ;    
}