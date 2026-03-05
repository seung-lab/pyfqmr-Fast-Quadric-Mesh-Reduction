#pragma once

/////////////////////////////////////////////
//
// Mesh Simplification
//
// (C) by Sven Forstmann in 2014
//
// License : MIT
// http://opensource.org/licenses/MIT
//
//https://github.com/sp4cerat/Fast-Quadric-Mesh-S;
// 5/2016: Chris Rorden created minimal version for OSX/Linux/Windows compile

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <span>
#include <string>
#include <vector>

#define loopi(start_l,end_l) for (auto i = start_l; i < end_l; ++i)
#define loopj(start_l,end_l) for (auto j = start_l; j < end_l; ++j)
#define loopk(start_l,end_l) for (auto k = start_l; k < end_l; ++k)

namespace Simplify {

constexpr size_t ZERO = 0;

struct vec3f {
	float x, y, z;

	inline vec3f( void ) {}

	inline vec3f( const float X, const float Y, const float Z )
	{ x = X; y = Y; z = Z; }

	inline vec3f operator + ( const vec3f& a ) const
	{ return vec3f( x + a.x, y + a.y, z + a.z ); }

	inline vec3f& operator += ( const vec3f& a )
	{ x += a.x; y += a.y; z += a.z; return *this; }

	inline vec3f operator * ( const float a ) const
	{ return vec3f( x * a, y * a, z * a ); }

	inline vec3f operator * ( const vec3f a ) const
	{ return vec3f( x * a.x, y * a.y, z * a.z ); }

	inline vec3f v3 () const
	{ return vec3f( x , y, z ); }

	inline vec3f& operator = ( const vec3f a )
	{ x=a.x;y=a.y;z=a.z;return *this; }

	inline vec3f operator / ( const vec3f a ) const
	{ return vec3f( x / a.x, y / a.y, z / a.z ); }

	inline vec3f operator - ( const vec3f& a ) const
	{ return vec3f( x - a.x, y - a.y, z - a.z ); }

	inline vec3f operator / ( const float a ) const
	{ return vec3f( x / a, y / a, z / a ); }

	inline float dot( const vec3f& a ) const
	{ return a.x*x + a.y*y + a.z*z; }

	inline vec3f cross( const vec3f& a , const vec3f& b )
	{
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
		return *this;
	}

	inline float angle( const vec3f& v )
	{
		vec3f a = v , b = *this;
		double dot = v.x*x + v.y*y + v.z*z;
		double len = a.length() * b.length();
		if(len==0)len=0.00001f;
		double input = dot  / len;
		if (input<-1) input=-1;
		if (input>1) input=1;
		return (float) acos ( input );
	}

	inline float angle2( const vec3f& v , const vec3f& w )
	{
		vec3f a = v, b = *this;
		double dot = a.x*b.x + a.y*b.y + a.z*b.z;
		double len = a.length() * b.length();
		if(len==0)len=1;

		vec3f plane; plane.cross( b,w );

		if ( plane.x * a.x + plane.y * a.y + plane.z * a.z > 0 ) {
			return (float) -acos ( dot  / len );
		}

		return (float) acos ( dot  / len );
	}

	inline vec3f rot_x( float a )
	{
		double yy = cos ( a ) * y + sin ( a ) * z;
		double zz = cos ( a ) * z - sin ( a ) * y;
		y = yy; z = zz;
		return *this;
	}
	inline vec3f rot_y( float a )
	{
		double xx = cos ( -a ) * x + sin ( -a ) * z;
		double zz = cos ( -a ) * z - sin ( -a ) * x;
		x = xx; z = zz;
		return *this;
	}
	inline void clamp( float min, float max )
	{
		if (x<min) x=min;
		if (y<min) y=min;
		if (z<min) z=min;
		if (x>max) x=max;
		if (y>max) y=max;
		if (z>max) z=max;
	}
	inline vec3f rot_z( double a )
	{
		double yy = cos ( a ) * y + sin ( a ) * x;
		double xx = cos ( a ) * x - sin ( a ) * y;
		y = yy; x = xx;
		return *this;
	}
	inline vec3f invert()
	{
		x=-x;y=-y;z=-z;return *this;
	}
	inline vec3f frac()
	{
		return vec3f(
			x-double(int(x)),
			y-double(int(y)),
			z-double(int(z))
			);
	}

	inline vec3f integer()
	{
	return vec3f(
		double(int(x)),
		double(int(y)),
		double(int(z))
		);
	}

	inline double length() const
	{
	return (double)sqrt(x*x + y*y + z*z);
	}

	inline vec3f normalize( double desired_length = 1 )
	{
		double square = sqrt(x*x + y*y + z*z);
		/*
		if (square <= 0.00001f )
		{
			x=1;y=0;z=0;
			return *this;
		}*/
		//double len = desired_length / square;
		x/=square;y/=square;z/=square;

		return *this;
	}
	static vec3f normalize( vec3f a );
};

vec3f barycentric(const vec3f &p, const vec3f &a, const vec3f &b, const vec3f &c) {
	vec3f v0 = b-a;
	vec3f v1 = c-a;
	vec3f v2 = p-a;
	double d00 = v0.dot(v0);
	double d01 = v0.dot(v1);
	double d11 = v1.dot(v1);
	double d20 = v2.dot(v0);
	double d21 = v2.dot(v1);
	double denom = d00*d11 - d01*d01;
	double v = (d11 * d20 - d01 * d21) / denom;
	double w = (d00 * d21 - d01 * d20) / denom;
	double u = 1.0 - v - w;
	return vec3f(u,v,w);
}

vec3f interpolate(
	const vec3f &p, 
	const vec3f &a, 
	const vec3f &b, 
	const vec3f &c, 
	const vec3f attrs[3]
) {
	vec3f bary = barycentric(p,a,b,c);
	vec3f out = vec3f(0,0,0);
	out = out + attrs[0] * bary.x;
	out = out + attrs[1] * bary.y;
	out = out + attrs[2] * bary.z;
	return out;
}

double min(double v1, double v2) {
	return fmin(v1,v2);
}

class SymmetricMatrix {
	public:
	float m[10];

	SymmetricMatrix(float c=0) { loopi(0,10) m[i] = c;  }

	SymmetricMatrix(
		float m11, float m12, float m13, float m14,
					float m22, float m23, float m24,
								float m33, float m34,
											float m44
	) {
		m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
								 m[4] = m22;  m[5] = m23;  m[6] = m24;
								 							m[7] = m33;  m[8] = m34;
																					 m[9] = m44;
	}

	// Make plane
	SymmetricMatrix(
		const float a, 
		const float b,
		const float c,
		const float d
	) {
		m[0] = a*a;  m[1] = a*b;  m[2] = a*c;  m[3] = a*d;
								 m[4] = b*b;  m[5] = b*c;  m[6] = b*d;
														  m[7] = c*c;  m[8] = c*d;
																					 m[9] = d*d;
	}

	double operator[](int c) const { return m[c]; }

	// Determinant
	double det(
		const int a11, const int a12, const int a13,
		const int a21, const int a22, const int a23,
		const int a31, const int a32, const int a33
	){
		return static_cast<double>(
			(m[a11]*m[a22]*m[a33]) + (m[a13]*m[a21]*m[a32]) + (m[a12]*m[a23]*m[a31])
			- (m[a13]*m[a22]*m[a31]) - (m[a11]*m[a23]*m[a32]) - (m[a12]*m[a21]*m[a33])
		);
	}

	const SymmetricMatrix operator+(const SymmetricMatrix& n) const {
		return SymmetricMatrix(
			m[0]+n[0],   m[1]+n[1],   m[2]+n[2],   m[3]+n[3],
			m[4]+n[4],   m[5]+n[5],   m[6]+n[6],
			m[7]+n[7],   m[8]+n[8],
			m[9]+n[9]
		);
	}

	SymmetricMatrix& operator+=(const SymmetricMatrix& n) {
		m[0]+=n[0];   m[1]+=n[1];   m[2]+=n[2];   m[3]+=n[3];
		m[4]+=n[4];   m[5]+=n[5];   m[6]+=n[6];   m[7]+=n[7];
		m[8]+=n[8];   m[9]+=n[9];
		return *this;
	}
};
///////////////////////////////////////////

// Global Variables & Structures
enum Attributes {
	NONE = 0,
	NORMAL = 2,
	TEXCOORD = 4,
	COLOR = 8
};

struct Triangle { 
	uint32_t v[3];
	float err[4];
	bool deleted;
	bool dirty;
	uint8_t attr;
	vec3f n;
	vec3f uvs[3];
	int16_t material; 
};
struct Vertex { 
	vec3f p;
	uint32_t tstart;
	uint32_t tcount;
	SymmetricMatrix q;
	bool border;
};

struct Mesh {
	std::vector<Triangle> triangles;
	std::vector<Vertex> vertices;

	std::string mtllib;
	std::vector<std::string> materials;

	void set(
		const std::span<float>& verts, 
		const std::span<uint32_t>& faces
	){
		vertices.clear();
		triangles.clear();

		uint64_t N_vertices = verts.size() / 3;
		uint64_t N_faces = faces.size() / 3;
		
		vertices.reserve(N_vertices * 3);
		triangles.reserve(N_faces * 3);

		for (uint64_t i = 0; i < N_vertices * 3; i += 3) {
			Vertex v;
			v.p.x = verts[i];
			v.p.y = verts[i+1];
			v.p.z = verts[i+2];
			vertices.push_back(v);
		}
		for (uint64_t i = 0; i < N_faces * 3; i += 3) {
			Triangle t;
			t.v[0] = faces[i];
			t.v[1] = faces[i+1];
			t.v[2] = faces[i+2];
			t.attr = 0;
			t.material = -1;
			triangles.push_back(t);
		}
	}

	std::vector<float> getVertices() const {
		std::vector<float> verts;
		const uint64_t N_vertices = vertices.size();
		verts.reserve(N_vertices * 3);
		loopi(0,N_vertices) {
			verts.push_back(vertices[i].p.x);
			verts.push_back(vertices[i].p.y);
			verts.push_back(vertices[i].p.z);
		}
		return verts;
	}

	std::vector<uint32_t> getFaces() const {
		std::vector<uint32_t> faces;
		const uint64_t N_faces = triangles.size();
		faces.reserve(N_faces * 3);
		loopi(0,N_faces) {
			faces.push_back(triangles[i].v[0]);
			faces.push_back(triangles[i].v[1]);
			faces.push_back(triangles[i].v[2]);
		}
		return faces;
	}

	std::vector<float> getNormals() {
		std::vector<float> normals;
		const uint64_t N_faces = triangles.size();
		normals.reserve(N_faces * 3);
		loopi(0, N_faces) {
			normals.push_back(triangles[i].n.x);
			normals.push_back(triangles[i].n.y);
			normals.push_back(triangles[i].n.z);
		}
		return normals;
	}
};

// Helper functions

double vertex_error(
	const SymmetricMatrix& q, 
	const double x, 
	const double y, 
	const double z
);
double calculate_error(
	const std::vector<Vertex>& vertices, 
	uint32_t id_v1, 
	uint32_t id_v2, 
	vec3f &p_result
);
bool flipped(
	const std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices, 
	const std::vector<uint32_t>& trefs,
	const std::vector<uint8_t>& vrefs,
	vec3f p,
	uint32_t i0,
	uint32_t i1,
	Vertex &v0,
	Vertex &v1,
	std::vector<int>& deleted
);
void update_uvs(
	std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices, 
	const std::vector<uint32_t>& trefs,
	const std::vector<uint8_t>& vrefs,
	uint32_t i0,
	const Vertex &v,
	const vec3f &p,
	std::vector<int> &deleted
);
void update_triangles(
	std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices,
	std::vector<uint32_t>& trefs,
	std::vector<uint8_t>& vrefs,
	uint32_t i0,
	Vertex& v,
	std::vector<int>& deleted,
	int64_t& deleted_triangles
);
void update_mesh(
	std::vector<Triangle>& triangles,
	std::vector<Vertex>& vertices,
	std::vector<uint32_t>& trefs,
	std::vector<uint8_t>& vrefs,
	int iteration
);
void compact_mesh(
	std::vector<Triangle>& triangles,
	std::vector<Vertex>& vertices
);

//
// Main simplification function
//
// target_count  : target nr. of triangles
// agressiveness : sharpness to increase the threshold.
//                 5..8 are good numbers
//                 more iterations yield higher quality
//
void simplify_mesh(
	Mesh& mesh,
	uint64_t target_count, 
	int update_rate = 5, 
	double agressiveness = 7,
	int max_iterations = 100,
	double alpha = 0.000000001,
	int K = 3,
	bool lossless = false,
	double threshold_lossless = 0.0001,
	bool preserve_border = false
) {
	std::vector<uint32_t> trefs;
	std::vector<uint8_t> vrefs;

	auto& triangles = mesh.triangles;
	auto& vertices = mesh.vertices;

	loopi(ZERO, triangles.size()) {
		triangles[i].deleted = 0;
	}

	// main iteration loop
	int64_t deleted_triangles = 0;
	std::vector<int> deleted0, deleted1;
	int64_t triangle_count = triangles.size();

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		if (triangle_count - deleted_triangles <= target_count) {
			break;
		}

		// update mesh once in a while
		if ((iteration % update_rate == 0) || lossless) {
			update_mesh(triangles, vertices, trefs, vrefs, iteration);
		}

		loopi(ZERO, triangles.size()) {
			triangles[i].dirty = 0;
		}

		//
		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		//
		double threshold = alpha * pow(double(iteration+K), agressiveness);
		if (lossless) {
			threshold = threshold_lossless;
		}

		// remove vertices & mark deleted triangles
		loopi(ZERO, triangles.size()) {
			Triangle &t = triangles[i];
			if (t.err[3] > threshold) continue;
			if (t.deleted) continue;
			if (t.dirty) continue;

			loopj(0,3) {
				if (t.err[j] < threshold) {
					auto i0 = t.v[ j     ];
					Vertex &v0 = vertices[i0];
					auto i1 = t.v[(j+1)%3];
					Vertex &v1 = vertices[i1];

					// Border check //Added preserve_border method from issue 14
					if (preserve_border) {
						if (v0.border || v1.border) continue; // should keep border vertices
					}
					else if (v0.border != v1.border) {
						continue; // base behaviour
					}

					// Compute vertex to collapse to
					vec3f p;
					calculate_error(vertices, i0,i1,p);
					deleted0.resize(v0.tcount); // normals temporarily
					deleted1.resize(v1.tcount); // normals temporarily
					// don't remove if flipped
					if (flipped(triangles,vertices,trefs,vrefs,p,i0,i1,v0,v1,deleted0)) {
						continue;
					}
					if (flipped(triangles,vertices,trefs,vrefs,p,i1,i0,v1,v0,deleted1)) {
						continue;
					}

					if ((t.attr & TEXCOORD) == TEXCOORD) {
						update_uvs(triangles,vertices,trefs,vrefs,i0,v0,p,deleted0);
						update_uvs(triangles,vertices,trefs,vrefs,i0,v1,p,deleted1);
					}

					// not flipped, so remove edge
					v0.p = p;
					v0.q = v1.q + v0.q;
					auto tstart = trefs.size();

					update_triangles(triangles,vertices, trefs, vrefs, i0, v0, deleted0, deleted_triangles);
					update_triangles(triangles,vertices, trefs, vrefs, i0, v1, deleted1, deleted_triangles);

					auto tcount = trefs.size() - tstart;

					if (tcount <= v0.tcount) {
						// save ram
						if (tcount) {
							memcpy(&trefs[v0.tstart], &trefs[tstart], tcount * sizeof(uint32_t));
							memcpy(&vrefs[v0.tstart], &vrefs[tstart], tcount * sizeof(uint8_t));
						}
					}
					else {
						// append
						v0.tstart = tstart;
					}

					v0.tcount=tcount;
					break;
				}
			}

			// done?
			if (lossless && (deleted_triangles <= 0)) {
				break;
			} 
			else if (!lossless && (triangle_count-deleted_triangles <= target_count)) {
				break;
			}

			if (lossless) deleted_triangles = 0;
		}
	}

	compact_mesh(triangles, vertices);
}

void simplify_mesh_lossless(
	Mesh& mesh,
	double epsilon = 1e-3,
	int max_iterations = 9999,
	bool preserve_border = false
) {
	std::vector<uint32_t> trefs;
	std::vector<uint8_t> vrefs;

	auto& triangles = mesh.triangles;
	auto& vertices = mesh.vertices;

	loopi(ZERO, triangles.size()) {
		triangles[i].deleted = 0;
	}
	// main iteration loop
	int64_t deleted_triangles = 0;
	std::vector<int> deleted0, deleted1;

	for (int iteration = 0; iteration < max_iterations; iteration++) {
		// update mesh constantly
		update_mesh(triangles, vertices, trefs, vrefs, iteration);
		// clear dirty flag
		loopi(ZERO, triangles.size()) triangles[i].dirty = 0;
		//
		// All triangles with edges below the threshold will be removed
		//
		// The following numbers works well for most models.
		// If it does not, try to adjust the 3 parameters
		//
		double threshold = epsilon;

		// remove vertices & mark deleted triangles
		loopi(ZERO, triangles.size()) {
			
			Triangle &t = triangles[i];
			if (t.err[3]>threshold) {
				continue;
			}
			else if (t.deleted) { 
				continue; 
			}
			else if (t.dirty) {
				continue;
			}

			loopj(0,3) {
				if (t.err[j] < threshold) {
					int i0=t.v[ j     ]; Vertex &v0 = vertices[i0];
					int i1=t.v[(j+1)%3]; Vertex &v1 = vertices[i1];

					// Border check //Added preserve_border method from issue 14 for lossless
					if (preserve_border) {
						if (v0.border || v1.border) continue; // should keep border vertices
					}
					else {
						if (v0.border != v1.border) continue; // base behaviour
					}

					// Compute vertex to collapse to
					vec3f p;
					calculate_error(vertices, i0, i1, p);

					deleted0.resize(v0.tcount); // normals temporarily
					deleted1.resize(v1.tcount); // normals temporarily

					// don't remove if flipped
					if (flipped(triangles,vertices,trefs,vrefs,p,i0,i1,v0,v1,deleted0)) continue;
					if (flipped(triangles,vertices,trefs,vrefs,p,i1,i0,v1,v0,deleted1)) continue;

					if ((t.attr & TEXCOORD) == TEXCOORD) {
						update_uvs(triangles,vertices,trefs,vrefs,i0,v0,p,deleted0);
						update_uvs(triangles,vertices,trefs,vrefs,i0,v1,p,deleted1);
					}

					// not flipped, so remove edge
					v0.p = p;
					v0.q = v1.q + v0.q;
					auto tstart = trefs.size();

					update_triangles(triangles,vertices,trefs,vrefs,i0,v0,deleted0,deleted_triangles);
					update_triangles(triangles,vertices,trefs,vrefs,i0,v1,deleted1,deleted_triangles);

					auto tcount = trefs.size() - tstart;

					if (tcount <= v0.tcount) {
						// save ram
						if (tcount) {
							memcpy(&trefs[v0.tstart], &trefs[tstart], tcount * sizeof(uint32_t));
							memcpy(&vrefs[v0.tstart], &vrefs[tstart], tcount * sizeof(uint8_t));
						}
					}
					else {
						v0.tstart = tstart; // append
					}

					v0.tcount = tcount;
					break;
				}
			}
		}
		
		if (deleted_triangles <= 0) {
			break;
		}
		
		deleted_triangles = 0;
	}
	compact_mesh(triangles, vertices);
}

// Check if a triangle flips when this edge is removed
bool flipped(
	const std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices,
	const std::vector<uint32_t>& trefs,
	const std::vector<uint8_t>& vrefs,
	vec3f p,
	uint32_t i0,
	uint32_t i1,
	Vertex& v0,
	Vertex& v1,
	std::vector<int>& deleted
) {
	loopk(0, v0.tcount) {
		const Triangle &t = triangles[trefs[v0.tstart+k]];
		if (t.deleted) continue;

		auto s = vrefs[v0.tstart+k];
		auto id1 = t.v[(s+1)%3];
		auto id2 = t.v[(s+2)%3];

		if (id1 == i1 || id2 == i1) { // delete ?
			deleted[k] = 1;
			continue;
		}
		vec3f d1 = vertices[id1].p-p;
		d1.normalize();
		vec3f d2 = vertices[id2].p-p;
		d2.normalize();
		
		if (fabs(d1.dot(d2)) > 0.999) {
			return true;
		}

		vec3f n;
		n.cross(d1,d2);
		n.normalize();
		deleted[k] = 0;
		if (n.dot(t.n) < 0.2) {
			return true;
		}
	}

	return false;
}

void update_uvs(
	std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices, 
	const std::vector<uint32_t>& trefs,
	const std::vector<uint8_t>& vrefs,
	uint32_t i0,
	const Vertex &v,
	const vec3f &p,
	std::vector<int> &deleted
) {
	loopk(0, v.tcount) {
		auto tid = trefs[v.tstart+k];
		auto tvertex = vrefs[v.tstart+k];

		Triangle &t = triangles[tid];
		if (t.deleted) {
			continue;
		}
		else if (deleted[k]) {
			continue;
		}

		vec3f p1 = vertices[t.v[0]].p;
		vec3f p2 = vertices[t.v[1]].p;
		vec3f p3 = vertices[t.v[2]].p;
		t.uvs[tvertex] = interpolate(p,p1,p2,p3,t.uvs);
	}
}

// Update triangle connections and edge error after a edge is collapsed
void update_triangles(
	std::vector<Triangle>& triangles,
	const std::vector<Vertex>& vertices, 
	std::vector<uint32_t>& trefs,
	std::vector<uint8_t>& vrefs,
	uint32_t i0,
	Vertex &v,
	std::vector<int> &deleted,
	int64_t &deleted_triangles
) {
	vec3f p;
	loopk(0, v.tcount) {
		auto tid = trefs[v.tstart+k];
		auto tvertex = vrefs[v.tstart+k];

		Triangle &t = triangles[tid];
		if (t.deleted) {
			continue;
		}
		else if (deleted[k]) {
			t.deleted = 1;
			deleted_triangles++;
			continue;
		}

		t.v[tvertex] = i0;
		t.dirty = 1;
		t.err[0] = calculate_error(vertices, t.v[0], t.v[1], p);
		t.err[1] = calculate_error(vertices, t.v[1], t.v[2], p);
		t.err[2] = calculate_error(vertices, t.v[2], t.v[0], p);
		t.err[3] = min(t.err[0], min(t.err[1],t.err[2]));
		
		trefs.push_back(tid);
		vrefs.push_back(tvertex);
	}
}

// compact triangles, compute edge error and build reference list
void update_mesh(
	std::vector<Triangle>& triangles,
	std::vector<Vertex>& vertices,
	std::vector<uint32_t>& trefs,
	std::vector<uint8_t>& vrefs,
	int iteration
) {
	if (iteration > 0) { // compact triangles
		int dst = 0;
		loopi(ZERO, triangles.size()) {
			if(!triangles[i].deleted) {
				triangles[dst++] = triangles[i];
			}
		}
		triangles.resize(dst);
	}

	// Init Reference ID list
	loopi(ZERO, vertices.size()) {
		vertices[i].tstart = 0;
		vertices[i].tcount = 0;
	}
	loopi(ZERO, triangles.size()) {
		Triangle &t = triangles[i];
		loopj(0,3) {
			vertices[t.v[j]].tcount++;
		}
	}

	uint32_t tstart = 0;
	loopi(ZERO, vertices.size()) {
		Vertex &v = vertices[i];
		v.tstart = tstart;
		tstart += v.tcount;
		v.tcount = 0;
	}

	// Write References
	trefs.resize(triangles.size()*3);
	vrefs.resize(triangles.size()*3);

	loopi(ZERO, triangles.size()) {
		Triangle &t = triangles[i];
		loopj(0,3) {
			Vertex &v = vertices[t.v[j]];
			trefs[v.tstart+v.tcount] = i;
			vrefs[v.tstart+v.tcount] = j;
			v.tcount++;
		}
	}

	// Identify boundary : vertices[].border=0,1
	if (iteration == 0) {
		std::vector<int> vcount,vids;

		loopi(ZERO, vertices.size()) {
			Vertex &v = vertices[i];
			vcount.clear();
			vids.clear();
			
			loopj(0, v.tcount) {
				uint32_t k = trefs[v.tstart+j];
				Triangle &t = triangles[k];
				loopk(0,3) {
					size_t ofs = 0;
					uint32_t id = t.v[k];
					while (ofs < vcount.size()) {
						if (vids[ofs] == id) {
							break;
						}
						ofs++;
					}
					if (ofs == vcount.size()) {
						vcount.push_back(1);
						vids.push_back(id);
					}
					else {
						vcount[ofs]++;
					}
				}
			}

			loopj(ZERO, vcount.size()) {
				vertices[vids[j]].border = (vcount[j] == 1);
			}
		}
	}

	//
	// Init Quadrics by Plane & Edge Errors
	//
	// required at the beginning ( iteration == 0 )
	// recomputing during the simplification is not required,
	// but mostly improves the result for closed meshes
	//
	if (iteration == 0) {
		loopi(ZERO, vertices.size()) {
			vertices[i].q = SymmetricMatrix(0.0);
		}

		loopi(ZERO, triangles.size()) {
			Triangle &t = triangles[i];
			vec3f n, p[3];
			loopj(0,3) {
				p[j] = vertices[t.v[j]].p;
			}
			n.cross(p[1]-p[0],p[2]-p[0]);
			n.normalize();
			t.n = n;
			loopj(0,3) {
				vertices[t.v[j]].q += SymmetricMatrix(n.x,n.y,n.z,-n.dot(p[0]));
			}
		}
		loopi(ZERO, triangles.size()) {
			// Calc Edge Error
			Triangle &t = triangles[i];vec3f p;
			loopj(0,3) t.err[j] = calculate_error(vertices, t.v[j], t.v[(j+1)%3], p);
			t.err[3] = min(t.err[0],min(t.err[1],t.err[2]));
		}
	}
}

// Finally compact mesh before exiting
void compact_mesh(
	std::vector<Triangle>& triangles,
	std::vector<Vertex>& vertices
) {
	int dst = 0;
	loopi(ZERO, vertices.size()) {
		vertices[i].tcount = 0;
	}
	loopi(ZERO, triangles.size()) {
		if(!triangles[i].deleted) {
			Triangle &t = triangles[i];
			triangles[dst++] = t;
			loopj(0,3) {
				vertices[t.v[j]].tcount = 1;
			}
		}
	}

	triangles.resize(dst);
	dst = 0;
	loopi(ZERO, vertices.size()) {
		if(vertices[i].tcount) {
			vertices[i].tstart=dst;
			vertices[dst].p=vertices[i].p;
			dst++;
		}
	}

	loopi(ZERO, triangles.size()) {
		Triangle &t = triangles[i];
		loopj(0,3)t.v[j] = vertices[t.v[j]].tstart;
	}
	vertices.resize(dst);
}

// Error between vertex and Quadric

double vertex_error(
	const SymmetricMatrix& q,
	const double x,
	const double y,
	const double z
) {
	return (
		q[0]*x*x 
		+ 2*q[1]*x*y 
		+ 2*q[2]*x*z 
		+ 2*q[3]*x 
		+ q[4]*y*y
		+ 2*q[5]*y*z
		+ 2*q[6]*y
		+ q[7]*z*z 
		+ 2*q[8]*z 
		+ q[9]
	);
}

// Error for one edge
double calculate_error(
	const std::vector<Vertex>& vertices,
	uint32_t id_v1, 
	uint32_t id_v2, 
	vec3f &p_result
) {
	// compute interpolated vertex

	SymmetricMatrix q = vertices[id_v1].q;
	q += vertices[id_v2].q;
	bool border = vertices[id_v1].border & vertices[id_v2].border;
	double error = 0;
	double det = q.det(0, 1, 2, 1, 4, 5, 2, 5, 7);
	if (det != 0 && !border) {
		// q_delta is invertible
		p_result.x = -1/det*(q.det(1, 2, 3, 4, 5, 6, 5, 7 , 8));  // vx = A41/det(q_delta)
		p_result.y =  1/det*(q.det(0, 2, 3, 1, 5, 6, 2, 7 , 8));  // vy = A42/det(q_delta)
		p_result.z = -1/det*(q.det(0, 1, 3, 1, 4, 6, 2, 5,  8));  // vz = A43/det(q_delta)
		error = vertex_error(q, p_result.x, p_result.y, p_result.z);
	}
	else {
		// det = 0 -> try to find best result
		vec3f p1 = vertices[id_v1].p;
		vec3f p2 = vertices[id_v2].p;
		vec3f p3 = (p1+p2)/2;
		double error1 = vertex_error(q, p1.x,p1.y,p1.z);
		double error2 = vertex_error(q, p2.x,p2.y,p2.z);
		double error3 = vertex_error(q, p3.x,p3.y,p3.z);
		error = min(error1, min(error2, error3));
		if (error1 == error) p_result = p1;
		if (error2 == error) p_result = p2;
		if (error3 == error) p_result = p3;
	}
	return (double)error;
}

char *trimwhitespace(char *str) {
	char *end;

	// Trim leading space
	while (isspace((unsigned char)*str)) {
		str++;
	}

	if (*str == 0) { // All spaces?
		return str;
	}

	// Trim trailing space
	end = str + strlen(str) - 1;
	while (end > str && isspace((unsigned char)*end)) {
		end--;
	}

	// Write new null terminator
	*(end+1) = 0;

	return str;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Not Useful for Trimesh wrapper, we rather need to directly assign vertices, faces colors etc. ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Mesh load_obj(const char* filename, bool process_uv = false) {
	Mesh mesh;

	mesh.vertices.clear();
	mesh.triangles.clear();
	
	if (filename == NULL || filename[0] == 0) {
		return mesh;
	}
	
	FILE* fn = fopen(filename, "rb");
	if (fn == NULL) {
		printf("File %s not found!\n", filename);
		return mesh;
	}
	
	fseek(fn, 0, SEEK_END);
	long file_size = ftell(fn);
	rewind(fn);
	
	mesh.vertices.reserve(file_size / 10);
	mesh.triangles.reserve(file_size / 10);
	
	char line[1000];
	int vertex_cnt = 0;
	int material = -1;
	std::map<std::string, int> material_map;
	std::vector<vec3f> uvs;
	std::vector<std::vector<int>> uvMap;
	
	while (fgets(line, sizeof(line), fn) != NULL) {
		char *p, *end;
		
		if (strncmp(line, "mtllib", 6) == 0) {
			mesh.mtllib = trimwhitespace(&line[7]);
		}
		else if (strncmp(line, "usemtl", 6) == 0) {
			std::string usemtl = trimwhitespace(&line[7]);
			if (material_map.find(usemtl) == material_map.end()) {
				material_map[usemtl] = mesh.materials.size();
				mesh.materials.push_back(usemtl);
			}
			material = material_map[usemtl];
		}
		else if (line[0] == 'v' && line[1] == 't' && line[2] == ' ') {
			vec3f uv;
			p = line + 3;
			uv.x = strtof(p, &end); if (end == p) continue; p = end;
			uv.y = strtof(p, &end); if (end == p) continue; p = end;
			uv.z = strtof(p, &end); // optional, ok if it fails
			if (end == p) uv.z = 0.0f;
			uvs.push_back(uv);
		}
		else if (line[0] == 'v' && line[1] == ' ') {
			Vertex v;
			p = line + 2;
			v.p.x = strtof(p, &end); if (end == p) continue; p = end;
			v.p.y = strtof(p, &end); if (end == p) continue; p = end;
			v.p.z = strtof(p, &end); if (end == p) continue;
			mesh.vertices.push_back(v);
		}
		else if (line[0] == 'f') {
			Triangle t;
			bool tri_ok = false;
			bool has_uv = false;
			int integers[9];
			
			if (sscanf(line,"f %d %d %d",
				&integers[0],&integers[1],&integers[2]) == 3) {
				tri_ok = true;
			}
			else if (sscanf(line,"f %d// %d// %d//",
				&integers[0],&integers[1],&integers[2]) == 3) {
				tri_ok = true;
			}
			else if (sscanf(line,"f %d//%d %d//%d %d//%d",
				&integers[0],&integers[3],
				&integers[1],&integers[4],
				&integers[2],&integers[5]) == 6) {
				tri_ok = true;
			}
			else if (sscanf(line,"f %d/%d/%d %d/%d/%d %d/%d/%d",
				&integers[0],&integers[6],&integers[3],
				&integers[1],&integers[7],&integers[4],
				&integers[2],&integers[8],&integers[5]) == 9) {
				tri_ok = true;
				has_uv = true;
			}
			else {
				printf("Unrecognized face: %s", line);
				fclose(fn);
				exit(0);
			}
			
			if (tri_ok) {
				t.v[0] = integers[0]-1-vertex_cnt;
				t.v[1] = integers[1]-1-vertex_cnt;
				t.v[2] = integers[2]-1-vertex_cnt;
				t.attr = 0;
				
				if (process_uv && has_uv) {
					std::vector<int> indices;
					indices.push_back(integers[6]-1-vertex_cnt);
					indices.push_back(integers[7]-1-vertex_cnt);
					indices.push_back(integers[8]-1-vertex_cnt);
					uvMap.push_back(indices);
					t.attr |= TEXCOORD;
				}

				t.material = material;
				mesh.triangles.push_back(t);
			}
		}
	}

	if (process_uv && uvs.size()) {
		loopi(ZERO, mesh.triangles.size()) {
			loopj(0,3) {
				mesh.triangles[i].uvs[j] = uvs[uvMap[i][j]];
			}
		}
	}

	fclose(fn);

	return mesh;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Not Useful for Trimesh wrapper, we rather need to directly return vertices, faces colors etc. ////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void write_obj(const Mesh& mesh, const char* filename) {
	auto& vertices = mesh.vertices;
	auto& triangles = mesh.triangles;
	auto& mtllib = mesh.mtllib;
	auto& materials = mesh.materials;

	FILE *file = fopen(filename, "w");
	int cur_material = -1;
	bool has_uv = (triangles.size() && (triangles[0].attr & TEXCOORD) == TEXCOORD);

	if (!file) {
		printf("write_obj: can't write data file \"%s\".\n", filename);
		exit(0);
	}
	if (!mtllib.empty()) {
		fprintf(file, "mtllib %s\n", mtllib.c_str());
	}
	loopi(ZERO, vertices.size()) {
		fprintf(file, "v %g %g %g\n", vertices[i].p.x,vertices[i].p.y,vertices[i].p.z); // more compact: remove trailing zeros
	}

	if (has_uv) {
		loopi(ZERO, triangles.size()) {
			if(!triangles[i].deleted) {
				fprintf(file, "vt %g %g\n", triangles[i].uvs[0].x, triangles[i].uvs[0].y);
				fprintf(file, "vt %g %g\n", triangles[i].uvs[1].x, triangles[i].uvs[1].y);
				fprintf(file, "vt %g %g\n", triangles[i].uvs[2].x, triangles[i].uvs[2].y);
			}
		}
	}

	int uv = 1;
	loopi(ZERO, triangles.size()) {
		if(!triangles[i].deleted) {
			if (triangles[i].material != cur_material) {
				cur_material = triangles[i].material;
				fprintf(file, "usemtl %s\n", materials[triangles[i].material].c_str());
			}
			if (has_uv) {
				fprintf(file, "f %d/%d %d/%d %d/%d\n", triangles[i].v[0]+1, uv, triangles[i].v[1]+1, uv+1, triangles[i].v[2]+1, uv+2);
				uv += 3;
			}
			else {
				fprintf(file, "f %d %d %d\n", triangles[i].v[0]+1, triangles[i].v[1]+1, triangles[i].v[2]+1);
			}
		}
	}

	fclose(file);
}

};
///////////////////////////////////////////

#undef loopi
#undef loopj
#undef loopk
