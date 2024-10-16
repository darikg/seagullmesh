from pathlib import Path

from seagullmesh import Mesh3
from seagullmesh import _seagullmesh as sgm

file = Path(__file__).parent / 'test' / 'assets' / 'armadillo.off'
m = Mesh3.from_file(str(file))

# m.interpolated_corrected_curvatures()
# print(m.vertex_data['principal_curvature'][m.vertices])

touched = m.vertex_data.add_property('touched', default=False)
vpm_mapper = sgm.meshing.VertexPointMapWrapper(m._mesh.points, touched.pmap)

x0 = sgm.meshing.UniformSizingField(0.5, m._mesh.points)
x1 = sgm.meshing.UniformSizingField(0.5, vpm_mapper)

# m.remesh()
# sizing = sgm.meshing.AdaptiveSizingField(tol=0.1, edge_len_min_max=(0.3, 0.5), faces=m.faces, mesh=m.mesh)