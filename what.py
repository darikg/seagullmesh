from pathlib import Path

from seagullmesh import Mesh3
from seagullmesh import _seagullmesh as sgm

file = Path(__file__).parent / 'test' / 'assets' / 'armadillo.off'
m = Mesh3.from_file(str(file))

# m.interpolated_corrected_curvatures()
# print(m.vertex_data['principal_curvature'][m.vertices])

print(m._mesh.points)
touched = m.vertex_data.add_property('touched', default=False)
print(sgm.meshing.VertexPointMapWrapper(m._mesh.points, touched.pmap))
# m.remesh()
# sizing = sgm.meshing.AdaptiveSizingField(tol=0.1, edge_len_min_max=(0.3, 0.5), faces=m.faces, mesh=m.mesh)