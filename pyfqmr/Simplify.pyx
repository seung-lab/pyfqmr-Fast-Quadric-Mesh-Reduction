# distutils: language = c++

cimport cython

from libcpp.vector cimport vector
from libcpp cimport bool

from time import time as _time

cimport numpy as np
import numpy as np

class _hidden_ref(object):
    """Hidden Python object to keep a reference to our numpy arrays
    """
    def __init__(self):
        self.faces = None
        self.verts = None

_REF = _hidden_ref() 

cdef extern from "Simplify.h":
    cdef cppclass vec3f:
       double x
       double y
       double z

cdef extern from "Simplify.h" namespace "Simplify" :
    void simplify_mesh( int target_count, int update_rate, double aggressiveness, 
                        bool verbose, int max_iterations,double alpha, int K, 
                        bool lossless, double threshold_lossless, bool preserve_border)
    void setMeshFromExt(vector[vector[double]] vertices, vector[vector[int]] faces)
    vector[vector[int]] getFaces()
    vector[vector[double]] getVertices()
    vector[vector[double]] getNormals()
    cdef cppclass Triangle:
        int v[3]
        int attr
        int material
    cdef cppclass Vertex:
        vec3f p
    vector[Triangle] triangles
    vector[Vertex] vertices


cdef class Simplify : 

    cdef vector[vector[int]] triangles_cpp
    cdef vector[vector[double]] vertices_cpp
    cdef vector[vector[double]] normals_cpp
    
    def __cinit__(self):
        pass

    @cython.boundscheck(False)
    def getMesh(self):
        """Gets the mesh from the simplify object once the simplification is done

        Returns
        -------
        verts : numpy.ndarray
            array of vertices of shape (n_vertices,3)
        faces : numpy.ndarray
            array of faces of shape (n_faces,3)
        norms : numpy.ndarray
            array of normals of shape (n_faces,3)
        """
        self.triangles_cpp = getFaces()
        self.vertices_cpp = getVertices()
        self.normals_cpp = getNormals()

        cdef size_t N_t = self.triangles_cpp.size()
        cdef size_t N_v = self.vertices_cpp.size()
        cdef size_t N_n = self.normals_cpp.size()
        cdef np.ndarray[int, ndim=2] faces = np.zeros((N_t, 3), dtype=np.int32)
        cdef np.ndarray[double, ndim=2] verts = np.zeros((N_v, 3), dtype=np.float64)
        cdef np.ndarray[double, ndim=2] norms = np.zeros((N_n, 3), dtype=np.float64)

        cdef size_t i = 0
        cdef size_t j = 0

        for i in range(N_t):
            for j in range(3):
                faces[i,j] = self.triangles_cpp[i][j]
        for i in range(N_v):
            for j in range(3):
                verts[i,j] = self.vertices_cpp[i][j]
        for i in range(N_n):
            for j in range(3):
                norms[i,j] = self.normals_cpp[i][j]

        return verts, faces, norms

    cpdef void setMesh(self, vertices_np, faces_np, face_colors=None):
        """Method to set the mesh of the simplifier object.
        
        Arguments
        ---------
        vertices : numpy.ndarray
            array of vertices of shape (n_vertices,3)
        faces : numpy.ndarray
            array of faces of shape (n_faces,3)
        face_colors : numpy.ndarray
            array of face_colors of shape (n_faces,3)
            this is not yet implemented
        """
        # Here we will need some checks, just to make sure the right objets are passed
        if vertices_np.dtype == np.float32:
            setVerticesNogil_float(vertices_np, vertices)
        elif vertices_np.dtype == np.float64:
            setVerticesNogil_double(vertices_np, vertices)
        else:
            setVerticesNogil_double(vertices_np.astype(dtype="float64", subok=False, copy=False), vertices)

        if faces_np.dtype == np.int32:
            setFacesNogil_int(faces_np, triangles)
        elif faces_np.dtype == np.uint32:
            setFacesNogil_uint(faces_np, triangles)
        else:
            setFacesNogil_int(faces_np.astype(dtype="int32", subok=False, copy=False), triangles)

    cpdef void simplify_mesh(self, int target_count = 100, int update_rate = 5, 
        double aggressiveness=7., max_iterations = 100, bool verbose=True,  
        bool lossless = False, double threshold_lossless=1e-3, double alpha = 1e-9, 
        int K = 3, bool preserve_border = True):
        """Simplify mesh

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
        t_start = _time()
        simplify_mesh(target_count, update_rate, aggressiveness, verbose, max_iterations, alpha, K,
                      lossless, threshold_lossless, preserve_border)
        t_end = _time()
        N_end = getFaces().size()


@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
cdef void setVerticesNogil_double(double[:,:] & vertices, vector[Vertex] & vector_vertices ) noexcept nogil:
    """nogil function for filling the vector of vertices, "vector_vertices",
    with the data found in the memory view of the array "vertices"
    """
    vector_vertices.resize(vertices.shape[0])

    cdef size_t i = 0
    for i in range(vertices.shape[0]):
        vector_vertices[i].p.x = vertices[i, 0];
        vector_vertices[i].p.y = vertices[i, 1];
        vector_vertices[i].p.z = vertices[i, 2];

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
cdef void setVerticesNogil_float(float[:,:] & vertices, vector[Vertex] & vector_vertices ) noexcept nogil:
    """nogil function for filling the vector of vertices, "vector_vertices",
    with the data found in the memory view of the array "vertices"
    """
    vector_vertices.resize(vertices.shape[0])

    cdef size_t i = 0
    for i in range(vertices.shape[0]):
        vector_vertices[i].p.x = vertices[i, 0];
        vector_vertices[i].p.y = vertices[i, 1];
        vector_vertices[i].p.z = vertices[i, 2];

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
cdef void setFacesNogil_int(int[:,:] & faces, vector[Triangle] & vector_faces ) noexcept nogil:
    """nogil function for filling the vector of faces, "vector_faces",
    with the data found in the memory view of the array "faces"
    """
    vector_faces.resize(faces.shape[0]);

    cdef size_t i = 0
    cdef size_t j = 0
    for i in range(faces.shape[0]):
        for j in range(3):
            vector_faces[i].v[j] = faces[i, j]

        vector_faces[i].attr = 0
        vector_faces[i].material = -1

@cython.boundscheck(False)
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)
cdef void setFacesNogil_uint(unsigned int[:,:] & faces, vector[Triangle] & vector_faces ) noexcept nogil:
    """nogil function for filling the vector of faces, "vector_faces",
    with the data found in the memory view of the array "faces"
    """
    vector_faces.resize(faces.shape[0]);

    cdef size_t i = 0
    cdef size_t j = 0
    for i in range(faces.shape[0]):
        for j in range(3):
            vector_faces[i].v[j] = faces[i, j]

        vector_faces[i].attr = 0
        vector_faces[i].material = -1




"""Example:

#We assume you have a numpy based mesh processing software
#Where you can get the vertices and faces of the mesh as numpy arrays.
#For example Trimesh or meshio
import pyfqmr
import trimesh as tr
bunny = tr.load_mesh('Stanford_Bunny_sample.stl')
#Simplify object
simp = pyfqmr.Simplify()
simp.setMesh(bunny.vertices, bunny.faces)
simp.simplify_mesh(target_count = 1000, aggressiveness=7, preserve_border=True, verbose=10)
vertices, faces, normals = simp.getMesh()
"""

"""Example2:
import trimesh as tr
import pyfqmr as fmr
mesh = tr.load_mesh('Stanford_Bunny_sample.stl')
simpl = fmr.Simplify()
verts, faces = mesh.vertices, mesh.faces
simpl.setMesh(verts, faces)
simpl.getMesh()
faces.shape
"""