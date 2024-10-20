import pytest
from pyvista import Cylinder, PolyData

from seagullmesh import Mesh3
from seagullmesh import _seagullmesh as sgm


def pv_border_edges(mesh: PolyData) -> PolyData:
    return mesh.extract_feature_edges(
        boundary_edges=True,
        non_manifold_edges=False,
        feature_edges=False,
        manifold_edges=False,
    )


@pytest.mark.skipif(not hasattr(sgm, 'border'), reason='border submodule not installed')
def test_label_border_vertices():
    pv_cyl = Cylinder(capping=False, resolution=10).triangulate()
    sm_cyl = Mesh3.from_pyvista(pv_cyl)

    sm_cyl.label_border_vertices('is_border')
    sm_n_verts = sm_cyl.vertex_data['is_border'][:].sum()
    pv_n_verts = pv_border_edges(pv_cyl).n_points

    assert sm_n_verts == pv_n_verts
