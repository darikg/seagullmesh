    typedef PMP::Adaptive_sizing_field<Mesh3, VertPoint>  AdaptiveSizingField;

    py::class_<AdaptiveSizingField>(sub, "AdaptiveSizingField", py::module_local())
        .def(py::init([](
                const double tol,
                const std::pair<double, double>& edge_len_min_max,
                const FaceRange& faces,
                Mesh3& mesh,
                double ball_radius
            ){
            auto params = PMP::parameters::.ball_radius(ball_radius);
            return AdaptiveSizingField(tol, edge_len_min_max, faces, mesh, params);
        })
    );
