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
	const int update_rate,
    const double aggressiveness,
    const int max_iterations,
    const bool verbose,
    const bool lossless,
    const double threshold_lossless,
    const double alpha,
    const int K,
    const bool preserve_border
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
		target_count, update_rate, aggressiveness,
		max_iterations, alpha, K, lossless, 
		threshold_lossless, preserve_border
	);

	if (verbose) {
		printf("finished simplifcation.\n");
	}

	const std::vector<float> decimated_vertices = mesh.getVertices();
	const std::vector<uint32_t> decimated_faces = mesh.getFaces();

	return py::make_tuple(
		py::array_t<float>(decimated_vertices.size(), decimated_vertices.data()),
		py::array_t<uint32_t>(decimated_faces.size(), decimated_faces.data())
	);
}

py::tuple simplify_from_obj(
	const std::string& filename,
	const bool process_uv,
	const uint64_t target_count, 
	const int update_rate,
    const double aggressiveness,
    const int max_iterations,
    const bool verbose,
    const bool lossless,
    const double threshold_lossless,
    const double alpha,
    const int K,
    const bool preserve_border
) {
	Simplify::Mesh mesh = Simplify::load_obj(filename.c_str(), process_uv);	

	if (verbose) {
		printf("loaded obj.\n");
	}

	Simplify::simplify_mesh(
		mesh,
		target_count, update_rate, aggressiveness,
		max_iterations, alpha, K, lossless, 
		threshold_lossless, preserve_border
	);

	if (verbose) {
		printf("finished simplifcation.\n");
	}

	const std::vector<float> decimated_vertices = mesh.getVertices();
	const std::vector<uint32_t> decimated_faces = mesh.getFaces();

	return py::make_tuple(
		py::array_t<float>(decimated_vertices.size(), decimated_vertices.data()),
		py::array_t<uint32_t>(decimated_faces.size(), decimated_faces.data())
	);
}

void simplify_from_obj_to_obj(
	const std::string& source,
	const std::string& dest,
	const bool process_uv,
	const int64_t target_count, 
	const int update_rate,
    const double aggressiveness,
    const int max_iterations,
    const bool verbose,
    const bool lossless,
    const double threshold_lossless,
    const double alpha,
    const int K,
    const bool preserve_border
) {
	Simplify::Mesh mesh = Simplify::load_obj(source.c_str(), process_uv);	

	if (verbose) {
		printf("loaded obj.\n");
	}

	Simplify::simplify_mesh(
		mesh,
		target_count, update_rate, aggressiveness,
		max_iterations, alpha, K, lossless, 
		threshold_lossless, preserve_border
	);

	if (verbose) {
		printf("finished simplifcation.\n");
		printf("writing obj.\n");
	}

	Simplify::write_obj(mesh, dest.c_str());

	if (verbose) {
		printf("done.\n");
	}
}

PYBIND11_MODULE(fastpyfqmr, m) {
	m.doc() = "Python Fast Quadratic Mesh Reduction. A quadrics based mesh simplification implementation."; 
	m.def("simplify_mesh", &simplify_mesh, "Apply quadrics based mesh simplification. Does not preserve topology.");
	m.def("simplify_from_obj", &simplify_from_obj, "Apply quadrics based mesh simplification. Does not preserve topology. Read directly from an obj.");
	m.def("simplify_from_obj_to_obj", &simplify_from_obj_to_obj, "Apply quadrics based mesh simplification. Does not preserve topology. Read directly from an obj and write directly to an obj.");
}
