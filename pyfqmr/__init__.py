import fastpyfqmr

import numpy as np
import numpy.typing as npt

def simplify_mesh(
	vertices:npt.NDArray[np.float32],
	faces:npt.NDArray[np.uint32],
	target_count:int = 100, 
	update_rate:int = 5,
	aggressiveness:float = 7.0,
	max_iterations:int = 100,
	verbose:bool = True,
	lossless:bool = False,
	threshold_lossless:float = 1e-3, 
	alpha:float = 1e-9,
	K:int = 3, 
	preserve_border:bool = True
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.uint32]]:
	"""
	Use quadrics based mesh decimation algorithm that does not preserve topology.

	Parameters
	----------
	target_count : int
	    Target number of triangles, not used if lossless is True
	update_rate : int
	    Number of iterations between each update.
	    If lossless flag is set to True, rate is 1
	aggressiveness : float
	    Parameter controlling the growth rate of the threshold at each
	    iteration when lossless is False.
	max_iterations : int
	    Maximal number of iterations
	verbose : bool
	    control verbosity
	lossless : bool
	    Use the lossless simplification method
	threshold_lossless : float
	    Maximal error after which a vertex is not deleted, only for
	    lossless method.
	alpha : float
	    Parameter for controlling the threshold growth
	K : int
	    Parameter for controlling the thresold growth
	preserve_border : Bool
	    Flag for preserving vertices on open border

	Note
	----
	threshold = alpha*pow( iteration + K, agressiveness)
	"""

	vertices = np.asarray(vertices, dtype=np.float32)
	vertices = vertices.reshape([ vertices.size, ], order="C")

	faces = np.asarray(faces, dtype=np.uint32)
	faces = faces.reshape([ faces.size, ], order="C")

	vertices, faces = fastpyfqmr.simplify_mesh(
		vertices, faces,
		target_count, update_rate, aggressiveness,
    	max_iterations, verbose, lossless,
    	threshold_lossless, alpha, K,
    	preserve_border,
	)

	vertices = vertices.reshape([ vertices.size // 3, 3 ], order="C")
	faces = faces.reshape([ faces.size // 3, 3 ], order="C")

	return vertices, faces

def simplify_mesh_from_obj(
	source:str,
	dest:Optiona[str] = None,
	target_count:int = 100, 
	update_rate:int = 5,
	aggressiveness:float = 7.0,
	max_iterations:int = 100,
	verbose:bool = True,
	lossless:bool = False,
	threshold_lossless:float = 1e-3, 
	alpha:float = 1e-9,
	K:int = 3, 
	preserve_border:bool = True
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.uint32]]:
	"""
	Use quadrics based mesh decimation algorithm that does not preserve topology.

	Parameters
	----------
	target_count : int
	    Target number of triangles, not used if lossless is True
	update_rate : int
	    Number of iterations between each update.
	    If lossless flag is set to True, rate is 1
	aggressiveness : float
	    Parameter controlling the growth rate of the threshold at each
	    iteration when lossless is False.
	max_iterations : int
	    Maximal number of iterations
	verbose : bool
	    control verbosity
	lossless : bool
	    Use the lossless simplification method
	threshold_lossless : float
	    Maximal error after which a vertex is not deleted, only for
	    lossless method.
	alpha : float
	    Parameter for controlling the threshold growth
	K : int
	    Parameter for controlling the thresold growth
	preserve_border : Bool
	    Flag for preserving vertices on open border

	Note
	----
	threshold = alpha*pow( iteration + K, agressiveness)
	"""

	vertices = np.asarray(vertices, dtype=np.float32)
	vertices = vertices.reshape([ vertices.size, ], order="C")

	faces = np.asarray(faces, dtype=np.uint32)
	faces = faces.reshape([ faces.size, ], order="C")

	vertices, faces = fastpyfqmr.simplify_mesh(
		vertices, faces,
		target_count, update_rate, aggressiveness,
    	max_iterations, verbose, lossless,
    	threshold_lossless, alpha, K,
    	preserve_border,
	)

	vertices = vertices.reshape([ vertices.size // 3, 3 ], order="C")
	faces = faces.reshape([ faces.size // 3, 3 ], order="C")

	return vertices, faces