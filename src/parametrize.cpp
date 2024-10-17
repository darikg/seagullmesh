#include "seagullmesh.hpp"

#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/LSCM_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/ARAP_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Two_vertices_parameterizer_3.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef Mesh3::Property_map<V, Point2>                                   UVMap;
typedef SMP::Two_vertices_parameterizer_3<Mesh3>                         TwoVertexBorderParameterizer;
typedef SMP::LSCM_parameterizer_3<Mesh3, TwoVertexBorderParameterizer>   LscmParameterizer;
typedef SMP::ARAP_parameterizer_3<Mesh3, TwoVertexBorderParameterizer>   ArapParameterizer;


void init_parametrize(py::module &m) {
    m.def_submodule("parametrize")
        .def("lscm", [](Mesh3& mesh, UVMap& uv) {
            H boundary_halfedge = PMP::longest_border(mesh).first;
            SMP::Error_code result = SMP::parameterize(mesh, LscmParameterizer(), boundary_halfedge, uv);
            return SMP::get_error_message(result);
        })
        .def("lscm", [](Mesh3& mesh, UVMap&uv, V v1, V v2) {
            H boundary_halfedge = PMP::longest_border(mesh).first;
            SMP::Error_code result = SMP::parameterize(
                mesh, LscmParameterizer(TwoVertexBorderParameterizer(v1, v2)), boundary_halfedge, uv);
            return SMP::get_error_message(result);
        })
        .def("arap", [](Mesh3& mesh, UVMap& uv) {
            H boundary_halfedge = PMP::longest_border(mesh).first;
            SMP::Error_code result = SMP::parameterize(mesh, ArapParameterizer(), boundary_halfedge, uv);
            return SMP::get_error_message(result);
        })
    ;
}