#pragma once

#include "AABBTree.h"

#include "OpenMeshCraft/Utils/SFINAE.h"

namespace OMC {

template <typename Traits>
class AABB_ProjectionTraversal
{
public:
	using NT      = typename Traits::NT;
	using PointT  = typename Traits::PointT;
	using PrimT   = typename Traits::PrimT;
	using BboxT   = typename Traits::BboxT;
	using SphereT = typename Traits::SphereT;

	using DoIntersect  = typename Traits::DoIntersect;
	using ProjectPoint = typename Traits::ProjectPoint;

	using PrimCPtr = const PrimT *;

public:
	AABB_ProjectionTraversal(const PointT &query, const PointT &hint,
	                         PrimCPtr hint_prim)
	  : m_query_sphere(query, (query - hint).sqrnorm())
	  , m_closest_point(hint)
	  , m_closest_prim(hint_prim)
	{
	}

	bool intersection(const PrimT &prim);

	bool do_inter(const BboxT &bbox) const;

	const PointT &closest_point() const { return m_closest_point; }
	PrimCPtr      closest_primitive() const { return m_closest_prim; }
	NT square_distance() const { return m_query_sphere.squared_radius(); }

private:
	SphereT      m_query_sphere;
	PointT       m_closest_point;
	PrimCPtr     m_closest_prim;
	DoIntersect  m_box_sphere_do_intersect;
	ProjectPoint m_project_point;
};

/// @todo Add an option to find first or all intersections.
template <typename Traits>
class AABB_BoxInterTraversal
{
public:
	using NT     = typename Traits::NT;
	using PointT = typename Traits::PointT;
	using PrimT  = typename Traits::PrimT;
	using QPrimT = typename Traits::QPrimT;
	using BboxT  = typename Traits::BboxT;

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

	using PrimCPtr  = const PrimT *;
	using PrimCPtrs = std::vector<PrimCPtr>;

public:
	AABB_BoxInterTraversal(const QPrimT &query)
	  : m_query(query)
	{
		m_box_of_query = m_calc_bbox(m_query);
	}

	bool intersection(const PrimT &prim);

	bool do_inter(const BboxT &bbox) const;

	const PrimCPtrs &result() const { return m_intersected_prims; }

private:
	CalcBbox    m_calc_bbox;
	DoIntersect m_do_intersect;

	QPrimT    m_query;
	BboxT     m_box_of_query;
	PrimCPtrs m_intersected_prims;
};

/// @todo Add an option to find first or all intersections.
template <typename Traits>
class AABB_PrimInterTraversal
{
public:
	using NT     = typename Traits::NT;
	using PointT = typename Traits::PointT;
	using PrimT  = typename Traits::PrimT;
	using QPrimT = typename Traits::QPrimT;
	using BboxT  = typename Traits::BboxT;

	using CalcBbox    = typename Traits::CalcBbox;
	using DoIntersect = typename Traits::DoIntersect;

	using PrimCPtr  = const PrimT *;
	using PrimCPtrs = std::vector<PrimCPtr>;

public:
	AABB_PrimInterTraversal(const QPrimT &query)
	  : m_query(query)
	{
		m_box_of_query = m_calc_bbox(m_query);
	}

	bool intersection(const PrimT &prim);

	bool do_inter(const BboxT &bbox) const;

	PrimCPtrs result() const { return m_intersected_prims; }

private:
	CalcBbox    m_calc_bbox;
	DoIntersect m_do_intersect;

	QPrimT    m_query;
	BboxT     m_box_of_query;
	PrimCPtrs m_intersected_prims;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
	#include "AABBTraversalTraits.inl"
#endif