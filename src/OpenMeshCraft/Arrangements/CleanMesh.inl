#pragma once

#include "CleanMesh.h"

namespace OMC {
template <typename Traits>
ArrCleanMesh<Traits>::ArrCleanMesh(const std::vector<NT>      &_in_coords,
                                   const std::vector<index_t> &_in_tris,
                                   const std::vector<size_t>  &_in_labels)
  : in_coords(_in_coords)
  , in_tris(_in_tris)
  , in_labels(_in_labels)
{
}

/// @brief Converts mesh identifiers (represented by indices) to labels
/// (represented by bitsets).
/// @pre No preconditions.
/// @post Calls `mergeDuplicatedVertices`.
template <typename Traits>
void ArrCleanMesh<Traits>::convertLabels()
{
	out_labels.resize(in_labels.size());
	Label mask;

	for (size_t i = 0; i < in_labels.size(); i++)
	{
		out_labels[i][in_labels[i]] = true;
		mask[in_labels[i]]          = true;
	}
	num_labels = mask.count();
}

/// @brief Merges duplicate vertices in `in_coords` and stores them in
/// `out_coords`. Fixes the indices in `in_tris` to produce `out_tris`.
/// @pre No required preconditions.
/// @post Calls `mergeDegeneratedAndDuplicatedTriangles`.
template <typename Traits>
void ArrCleanMesh<Traits>::mergeDuplicatedVertices()
{
	bool parallel = in_coords.size() / 3 > 10000;

	// reserve memory for output
	out_coords.reserve(in_coords.size());
	// define a temporary vec3 type for sorting vertices.
	using vec3                      = std::array<NT, 3>;
	vec3                *in_vecs    = (vec3 *)in_coords.data();
	// record number of vertices
	size_t               origin_num = in_coords.size() / 3;
	size_t               unique_num = 0;
	// define a temporary vector of vertex indices for sort.
	std::vector<index_t> sorted(origin_num);
	std::iota(sorted.begin(), sorted.end(), 0);
	// sort vertices
	if (parallel)
		tbb::parallel_sort(sorted.begin(), sorted.end(), [in_vecs](auto a, auto b)
		                   { return in_vecs[a] < in_vecs[b]; });
	else
		std::sort(sorted.begin(), sorted.end(),
		          [in_vecs](auto a, auto b) { return in_vecs[a] < in_vecs[b]; });
	// use a lookup table to remove duplicate vertices and update their indices.
	std::vector<index_t> lookup(origin_num);
	for (size_t idx = 0; idx < origin_num; idx++)
	{
		if (idx == 0 || in_vecs[sorted[idx]] != in_vecs[sorted[idx - 1]])
		{
			const vec3 &v = in_vecs[sorted[idx]];
			// save the unique vertex into out_coords
			out_coords.push_back(v[0]);
			out_coords.push_back(v[1]);
			out_coords.push_back(v[2]);
			unique_num += 1;
		}
		lookup[sorted[idx]] = unique_num - 1;
	}
	// update vertex indices in triangles.
	out_tris.resize(in_tris.size());
	if (parallel)
		std::transform(std::execution::par_unseq, in_tris.begin(), in_tris.end(),
		               out_tris.begin(),
		               [&lookup](index_t idx) { return lookup[idx]; });
	else
		std::transform(std::execution::seq, in_tris.begin(), in_tris.end(),
		               out_tris.begin(),
		               [&lookup](index_t idx) { return lookup[idx]; });
}

/// @brief Removes degenerate and duplicate triangles from `out_tris` and
/// merges the labels of duplicate triangles into a unique label in
/// `out_labels`. As a result, multiple bits in `out_labels` (std::bitset) may
/// be set.
/// @pre Calls `removeDuplicatedVertices`.
/// @post Calls `mergeDegeneratedAndDuplicatedTriangles`.
template <typename Traits>
void ArrCleanMesh<Traits>::removeDegenerateAndDuplicatedTriangles()
{
	using vec3            = std::array<NT, 3>;
	using vec3i           = std::array<index_t, 3>;
	size_t num_orig_verts = out_coords.size() / 3;
	size_t num_orig_tris  = out_tris.size() / 3;
	vec3  *ptr_verts      = (vec3 *)out_coords.data();
	vec3i *ptr_tris       = (vec3i *)out_tris.data();
	size_t tri_off        = 0;

	// compute collinear
	std::vector<uint8_t> collinear_res(num_orig_tris, false);
	tbb::parallel_for((size_t)0, num_orig_tris,
	                  [this, ptr_verts, ptr_tris, &collinear_res](index_t t_id)
	                  {
		                  const vec3i &t      = ptr_tris[t_id];
		                  collinear_res[t_id] = CollinearPoints3D()(
		                    ptr_verts[t[0]].data(), ptr_verts[t[1]].data(),
		                    ptr_verts[t[2]].data());
	                  });

	bool parallel = num_orig_tris > 10000;
	if (parallel)
	{
		// map: edge (two larger tri vertices) -> triangle_id
		using TriMap = phmap::flat_hash_map<std::pair<index_t, index_t>, index_t>;
		// the smallest tri vertex -> map
		std::vector<tbb::spin_mutex> tris_map_mutex(num_orig_verts);
		std::vector<TriMap>          tris_map(num_orig_verts);

		// unique index or compact index
		std::vector<index_t> tris_idx(num_orig_tris);

		bool exist_removed_tri = false;

		auto check_tri_is_unique = [this, &collinear_res, &tris_map,
		                            &tris_map_mutex, &tris_idx,
		                            &exist_removed_tri](index_t t_id)
		{
			if (collinear_res[t_id])
			{
				tris_idx[t_id]    = InvalidIndex; // means removed
				exist_removed_tri = true;
				return;
			}

			index_t v0_id = out_tris[(3 * t_id)];
			index_t v1_id = out_tris[(3 * t_id) + 1];
			index_t v2_id = out_tris[(3 * t_id) + 2];

			vec3i tri = {v0_id, v1_id, v2_id};
			std::sort(tri.begin(), tri.end());

			std::pair<index_t, index_t> edge = {tri[1], tri[2]};

			{ // critical section
				std::lock_guard<tbb::spin_mutex> lock(tris_map_mutex[tri[0]]);

				auto ins = tris_map[tri[0]].insert({edge, t_id});
				// if we meet this triangle for the first time, it is unique and
				// ins.first->second is same as t_id. otherwise, it is duplicate and
				// ins.first->second is the index of existed triangle.
				if (t_id < ins.first->second)
				{ // always save unique triangle with lower index
					tris_idx[ins.first->second] = t_id;
					tris_idx[t_id]              = t_id;
					tris_map[tri[0]][edge]      = t_id;
				}
				else
					tris_idx[t_id] = ins.first->second;
				if (!ins.second)
					exist_removed_tri = true;
			}
		};

		// parallel build map from triangle to unique tri
		tbb::parallel_for((size_t)0, num_orig_tris, check_tri_is_unique);

		if (!exist_removed_tri)
		{
			// no triangles were removed; exiting the function.
			return;
		}
		// clear to save memory
		tris_map       = std::vector<TriMap>();
		tris_map_mutex = std::vector<tbb::spin_mutex>();

		// loop as before by a simpler way
		std::vector<index_t> compact_tris_idx(num_orig_tris);

		// the first traversal to map unique triangle to compact index
		for (index_t t_id = 0; t_id < num_orig_tris; t_id++)
		{
			if (!is_valid_idx(tris_idx[t_id]))
				continue;

			index_t v0_id = out_tris[3 * t_id];
			index_t v1_id = out_tris[3 * t_id + 1];
			index_t v2_id = out_tris[3 * t_id + 2];
			Label   label = out_labels[t_id];

			if (tris_idx[t_id] == t_id)
			{
				// triangle is unique, save it and its compact id.
				out_tris[tri_off * 3]     = v0_id;
				out_tris[tri_off * 3 + 1] = v1_id;
				out_tris[tri_off * 3 + 2] = v2_id;
				out_labels[tri_off]       = label;

				// update unique index to compact index
				compact_tris_idx[t_id] = tri_off;
				tri_off += 1;
			}
			else
			{
				// triangle is duplicate, save info abount duplication
				index_t unique_idx = t_id;
				while (unique_idx != tris_idx[unique_idx])
					unique_idx = tris_idx[unique_idx];
				index_t compact_idx = compact_tris_idx[unique_idx];

				out_labels[compact_idx] |= label;
				size_t mesh_label = LabelToIdx(label);

				index_t curr_tri_verts[] = {v0_id, v1_id, v2_id};
				index_t uniq_tri_verts[] = {out_tris[compact_idx * 3],
				                            out_tris[compact_idx * 3 + 1],
				                            out_tris[compact_idx * 3 + 2]};
				bool    w = consistentWinding(curr_tri_verts, uniq_tri_verts);
				dupl_triangles.push_back({/*compact triangle id*/ compact_idx,
				                          /*label of the actual triangle*/ mesh_label,
				                          /*winding with respect to the triangle stored
				                             in mesh (true -> same, false -> opposite)*/
				                          w});
			}
		}
	}
	else
	{
		// loop as before by use simpler way
		// map: tri_vertices -> tri_off
		phmap::flat_hash_map<vec3i, size_t, hash<vec3i>> tris_map;

		for (size_t t_id = 0; t_id < num_orig_tris; ++t_id)
		{
			if (collinear_res[t_id])
				continue;
			index_t v0_id = out_tris[(3 * t_id)];
			index_t v1_id = out_tris[(3 * t_id) + 1];
			index_t v2_id = out_tris[(3 * t_id) + 2];
			Label   label = out_labels[t_id];

			vec3i tri = {v0_id, v1_id, v2_id};
			std::sort(tri.begin(), tri.end());

			auto ins = tris_map.insert({tri, tri_off});

			if (ins.second) // first time for tri v0, v1, v2
			{
				out_tris[tri_off * 3]     = v0_id;
				out_tris[tri_off * 3 + 1] = v1_id;
				out_tris[tri_off * 3 + 2] = v2_id;
				out_labels[tri_off]       = label;
				tri_off += 1;
			}
			else // triangle already present -> save info about duplicates
			{
				size_t orig_tri_off = ins.first->second;
				out_labels[orig_tri_off] |= label; // label for duplicates

				size_t mesh_label = LabelToIdx(label);
				OMC_EXPENSIVE_ASSERT(mesh_label >= 0, "invalid label id");

				index_t curr_tri_verts[] = {v0_id, v1_id, v2_id};
				index_t orig_tri_verts[] = {out_tris[orig_tri_off * 3],
				                            out_tris[orig_tri_off * 3 + 1],
				                            out_tris[orig_tri_off * 3 + 2]};

				bool w = consistentWinding(curr_tri_verts, orig_tri_verts);

				dupl_triangles.push_back({/*original triangle id*/ orig_tri_off,
				                          /*label of the actual triangle*/ mesh_label,
				                          /*winding with respect to the triangle stored
				                             in mesh (true -> same, false -> opposite)*/
				                          w});
			}
		}
	}

	out_tris.resize(tri_off * 3);
	out_labels.resize(tri_off);
}

/// @brief After removing degenerate triangles, there are possibly isolated
/// vertices, so remove them.
/// @pre Calls `removeDegenerateAndDuplicatedTriangles`
/// @post No required postconditions.
template <typename Traits>
void ArrCleanMesh<Traits>::removeIsolatedVertices()
{
	std::vector<size_t> vert_mark(out_coords.size() / 3, /*true*/ 1);
	// check if vertices are isolated.
	for (index_t vi : out_tris)
		vert_mark[vi] = /*false*/ 0; // set it to false if it is not isolated.
		                             // remove isolated vertices and update index
	size_t num_orig_verts = out_coords.size() / 3;
	size_t num_conn_verts = 0;
	for (index_t vj = 0; vj < num_orig_verts; vj++)
	{
		if (!vert_mark[vj])
		{
			// if vertex `vj` is not isolated, store its new index `vi`
			index_t vi = num_conn_verts;
			num_conn_verts += 1;

			vert_mark[vj]          = vi;
			out_coords[vi * 3]     = out_coords[vj * 3];
			out_coords[vi * 3 + 1] = out_coords[vj * 3 + 1];
			out_coords[vi * 3 + 2] = out_coords[vj * 3 + 2];
		}
	}

	if (num_orig_verts == num_conn_verts) // no isolated vertex is found
		return;

	// update vertex index in out_tris
	for (index_t &vi : out_tris)
		vi = vert_mark[vi];
}

/// @brief t0 -> vertex ids of triangle t0, t1 -> vertex ids of triangle t1
template <typename Traits>
bool ArrCleanMesh<Traits>::consistentWinding(const index_t *t0,
                                             const index_t *t1)
{
	int j = 0;
	while (j < 3 && t0[0] != t1[j])
		j++;
	OMC_EXPENSIVE_ASSERT(j < 3, "not same triangle");
	OMC_EXPENSIVE_ASSERT((t0[1] == t1[(j + 1) % 3] && t0[2] == t1[(j + 2) % 3]) ||
	                       (t0[1] == t1[(j + 2) % 3] && t0[2] == t1[(j + 1) % 3]),
	                     "not same triangle");
	return t0[1] == t1[(j + 1) % 3] && t0[2] == t1[(j + 2) % 3];
}

} // namespace OMC