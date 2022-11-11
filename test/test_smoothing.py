import pytest
from pyvista import Sphere

from seagullmesh import Mesh3


@pytest.fixture()
def sphere():
    return Mesh3.from_pyvista(Sphere().triangulate())


@pytest.mark.parametrize(('angle', 'area'), [(True, False), (False, True), (True, True)])
def test_angle_and_area_smoothing(sphere, area, angle):
    sphere.smooth_angle_and_area(
        sphere.faces,
        n_iter=1,
        use_area_smoothing=area,
        use_angle_smoothing=angle,
    )


def test_tangential_relaxation(sphere):
    sphere.tangential_relaxation(sphere.vertices, n_iter=1)


def test_smooth_shape(sphere):
    sphere.smooth_shape(sphere.faces, n_iter=1, time=.1)
