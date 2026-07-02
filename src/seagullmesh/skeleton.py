from __future__ import annotations

from functools import cached_property
from typing import Dict, List

import numpy as np

from seagullmesh import Mesh3


class Skeleton:
    """Wrapper around the C++ sgm.skeletonization.Skeleton class

    (Which is itself a boost adjacency_list)
    """
    def __init__(self, mesh: Mesh3, skeleton):
        self.mesh = mesh
        self._skeleton = skeleton

    @cached_property
    def points(self) -> np.ndarray:
        return self._skeleton.points

    @cached_property
    def edges(self) -> np.ndarray:
        return self._skeleton.edges

    @cached_property
    def vertex_map(self) -> Dict[int, List[Vertex]]:
        return self._skeleton.vertex_map

    @cached_property
    def radii(self) -> np.ndarray:
        return self._skeleton.compute_radii(self.mesh.mesh)

    def to_pyvista(self) -> pv.PolyData:
        import pyvista as pv
        sk_mesh = pv.PolyData()
        sk_mesh.points = self.points
        sk_mesh.lines = pv.CellArray.from_regular_cells(self.edges)
        sk_mesh.point_data['min_radius'] = self.radii[:, 0]
        sk_mesh.point_data['max_radius'] = self.radii[:, 1]
        return sk_mesh
