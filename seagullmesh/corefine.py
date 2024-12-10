from dataclasses import dataclass
from typing_extensions import Self

from seagullmesh import Mesh3, PropertyMap, Vertex, Edge, sgm
from seagullmesh._seagullmesh import corefine


@dataclass
class InputSpec:
    mesh: Mesh3
    vertex_index: PropertyMap[Vertex, int] | None = None
    edge_constrained: PropertyMap[Edge, bool] = None
    face_index: PropertyMap[Vertex, int] | None = None


class Corefiner:
    def __init__(
            self,
            mesh0: Mesh3,
            mesh1: Mesh3,
            inplace=False,
    ):
        self.specs = InputSpec(mesh=mesh0), InputSpec(mesh=mesh1)
        self.output = mesh0 if inplace else Mesh3()

    def vertex_index(self, name: str) -> Self:
        for s in self.specs:
            s.vertex_index = s.mesh.vertex_data.get_or_create_property(name, default=-1)  # TODO dtype signed_int
        return self

    def edge_constrained(self, name: str) -> Self:
        for s in self.specs:
            s.edge_constrained = s.mesh.edge_data.get_or_create_property(name, default=False)
        return self

    def face_index(self, name: str) -> Self:
        for s in self.specs:
            s.face_index = s.mesh.face_data.get_or_create_property(name, default=-1)
        return self

    def visitor(self) -> corefine.CorefinementVertexTracker | corefine.CorefinementVertexFaceTracker | None:
        if any(s.face_index for s in self.specs):
            # tracker = sgm.corefine.CorefinementVertexFaceTracker(
            #     mesh1.mesh, mesh2.mesh, vert_idx1.pmap, vert_idx2.pmap, face_idx1.pmap, face_idx2.pmap)
            pass # return face tracker

        if any(s.vertex_index for s in self.specs):
            # tracker = sgm.corefine.CorefinementVertexTracker(mesh1.mesh, mesh2.mesh, vert_idx1.pmap, vert_idx2.pmap)
            pass # return vertex tracker

        return None

    def corefine(self, other: Mesh3) -> None:
        """Corefines the two meshes in place"""
        sgm.corefine.corefine(self._mesh, other._mesh)


    def union(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean union"""

        sgm.corefine.union(self._mesh, other._mesh, out._mesh)
        return out


    def difference(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean difference"""
        out = self if inplace else Mesh3(_Mesh3())
        sgm.corefine.difference(self._mesh, other._mesh, out._mesh)
        return out


    def intersection(self, other: Mesh3, inplace=False) -> Mesh3:
        """Corefines the two meshes and returns their boolean intersection"""
        out = self if inplace else Mesh3(_Mesh3())
        sgm.corefine.intersection(self._mesh, other._mesh, out._mesh)
        return out
