#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <span>

#include "Simplify.h"

namespace py = pybind11;

py::tuple simplify_mesh(
	const py::array_t<float>& vertices,
	const py::array_t<uint32_t>& faces,
	const uint64_t target_count, 
	const int update_rate = 5,
    const double aggressiveness = 7.0, 
    const int max_iterations = 100, 
    const bool verbose = true,
    const bool lossless = false, 
    const double threshold_lossless = 1e-3,
    const double alpha = 1e-9,
    const int K = 3, 
    const bool preserve_border = true
) {
	Simplify::Mesh mesh;

	auto vert_buf = vertices.request();
	auto face_buf = faces.request();

	mesh.set(
		std::span<float>((float*)vert_buf.ptr, vert_buf.size),
		std::span<uint32_t>((uint32_t*)face_buf.ptr, face_buf.size)
	);

	Simplify::simplify_mesh(
		mesh,
		target_count, update_rate, agressiveness,
		max_iterations, alpha, K, lossless, 
		threshold_lossless, preserve_border
	);

	const std::vector<float> decimated_vertices = mesh.getVertices();
	const std::vector<uint32_t> decimated_faces = mesh.getFaces();

	return py::make_tuple(
		py::array_t<float>(decimated_vertices.size(), decimated_vertices.data()),
		py::array_t<uint32_t>(decimated_faces.size(), decimated_faces.data())
	);
}

PYBIND11_MODULE(fastpyfqmr, m) {
	m.doc() = "Python Fast Quadratic Mesh Reduction. A quadrics based mesh simplification implementation."; 
	m.def("simplify_mesh", &simplify_mesh, "Apply quadrics based mesh simplification. Does not preserve topology.");
}
