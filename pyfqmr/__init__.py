import fastpyfqmr

import numpy as np

def simplify(vertices, faces):
	
	vertices = np.asarray(vertices, dtype=np.float32)
	vertices = vertices.reshape([ vertices.size, ], order="C")

	faces = np.asarray(faces, dtype=np.uint32)
	faces = faces.reshape([ faces.size, ], order="C")

	vertices, faces = fastpyfqmr.simplify_mesh(vertices, faces)

	vertices = vertices.reshape([ vertices.size // 3, 3 ], order="C")
	faces = faces.reshape([ faces.size // 3, 3 ], order="C")

	return vertices, faces