#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>

typedef PMP::Principal_curvatures_and_directions<Kernel>    PrincipalCurvDir;

typedef Mesh3::Property_map<V, double>                      VertDouble;
typedef Mesh3::Property_map<V, PrincipalCurvDir>            VertPrincipalCurvDir;

void init_meshing(py::module &m) {
    m.def_submodule("meshing")
        .def("interpolated_corrected_curvatures", [](
            const Mesh3& mesh,
            const double ball_radius,
            VertDouble& mean_curv_map,
            VertDouble& gauss_curv_map,
            VertPrincipalCurvDir princ_curv_dir_map,
        ) {
            auto params = PMP::parameters::
                vertex_mean_curvature_map(mean_curv_map)
                .vertex_Gaussian_curvature_map(mean_curv_map)
                .vertex_principal_curvatures_and_directions_map(princ_curv_dir_map)
                .ball_radius(ball_radius)
            ;
            PMP::interpolated_corrected_curvatures(mesh, params);
        })

    ;
}
