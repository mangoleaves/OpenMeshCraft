#pragma once

#include "KdTree.h"

namespace OMC {

template <typename KdTraits>
template <typename Points, typename Attrs>
void KdTree<KdTraits>::insert(const Points &points, const Attrs &attrs)
{
	static_assert(
	  std::is_same_v<remove_cvref_t<decltype(*points.begin())>, PointT>,
	  "Points' types are different.");

	static_assert(
	  std::is_same_v<remove_cvref_t<decltype(*attrs.begin())>, PointAttrT>,
	  "Attributes' types are different.");

	OMC_THROW_INVALID_ARGUMENT_IF(points.size() != attrs.size(),
	                           "size of points and attrs are different.");

	pts.reserve(points.size());

	auto pnt_first  = points.begin();
	auto attr_first = attrs.begin();
	for (; pnt_first != points.end(); pnt_first++, attr_first++)
		pts.emplace_back(*pnt_first, *attr_first);
}

template <typename KdTraits>
void KdTree<KdTraits>::build()
{
	data.reserve(pts.size());
	for (size_t i = 0; i < pts.size(); i++)
	{
		data.emplace_back(&pts[i]);
	}

	PointContainer container(data.begin(), data.end());
	bbox = container.bbox();
	if (container.size() <= bucket_size)
	{
		tree_root = create_leaf_node(container);
	}
	else
	{
		tree_root = new_internal_node();
		create_internal_node(tree_root, container);
	}

	// Sort points for better memory cache.
	KdPoints pts_tmp;
	pts_tmp.reserve(pts.size());
	for (size_t i = 0; i < pts.size(); i++)
	{
		pts_tmp.emplace_back(*data[i]);
	}
	for (size_t i = 0; i < leaf_nodes.size(); i++)
	{
		std::ptrdiff_t tmp = leaf_nodes[i].data - pts.begin();
		leaf_nodes[i].data = pts_tmp.begin() + tmp;
	}
	pts.swap(pts_tmp);

	data.clear();
}

template <typename KdTraits>
[[nodiscard]] auto KdTree<KdTraits>::search_nearest_point(const PointT &query)
  -> std::pair<PointT, PointAttrT>
{
	OrthogonalNearestSeach<KdTraits> ons(tree_root, bbox, query);
	KdPointsIter                     kp_it = ons.closest_point_iter();
	return std::pair<PointT, PointAttrT>(kp_it->primitive(), kp_it->attribute());
}

template <typename KdTraits>
auto KdTree<KdTraits>::create_leaf_node(PointContainer &container) -> NodePtr
{
	leaf_nodes.emplace_back(container.size());

	std::ptrdiff_t tmp = container.begin() - data.begin();

	leaf_nodes.back().data = pts.begin() + tmp;

	return &(leaf_nodes.back());
}

template <typename KdTraits>
auto KdTree<KdTraits>::new_internal_node() -> NodePtr
{
	internal_nodes.emplace_back();
	return &(internal_nodes.back());
}

template <typename KdTraits>
void KdTree<KdTraits>::create_internal_node(NodePtr         n,
                                            PointContainer &container)
{
	InternalPtr nh = static_cast<InternalPtr>(n);

	size_t         sep_dim;
	NT             sep_val;
	PointContainer container_low;

	split(sep_dim, sep_val, container, container_low);
	nh->set_separator(sep_dim, sep_val);

	handle_extended_node(nh, container, container_low);

	if (container_low.size() > bucket_size)
	{
		nh->m_lower_ch = new_internal_node();
		create_internal_node(nh->m_lower_ch, container_low);
	}
	else
	{
		nh->m_lower_ch = create_leaf_node(container_low);
	}

	if (container.size() > bucket_size)
	{
		nh->m_upper_ch = new_internal_node();
		create_internal_node(nh->m_upper_ch, container);
	}
	else
	{
		nh->m_upper_ch = create_leaf_node(container);
	}
}

template <typename KdTraits>
void KdTree<KdTraits>::split(size_t &sep_dim, NT &sep_val,
                             PointContainer &container_origin,
                             PointContainer &container_low)
{
	sep_dim = container_origin.bbox().longest_axis();

	// Fix degenerated cases
	if (container_origin.tbox().min_coord(sep_dim) !=
	    container_origin.tbox().max_coord(sep_dim))
	{
		sep_val = (container_origin.bbox().max_coord(sep_dim) +
		           container_origin.bbox().min_coord(sep_dim)) /
		          2.0;
	}
	else
	{
		sep_dim = container_origin.tbox().longest_axis();
		sep_val = (container_origin.tbox().max_coord(sep_dim) +
		           container_origin.tbox().min_coord(sep_dim)) /
		          2.0;
	}

	const NT &max_span_upper = container_origin.tbox().max_coord(sep_dim);
	const NT &max_span_lower = container_origin.tbox().min_coord(sep_dim);
	if (max_span_upper <= sep_val)
	{
		sep_val = max_span_upper;
	}
	if (max_span_lower >= sep_val)
	{
		sep_val = max_span_lower;
	}

	container_origin.split(container_low, sep_dim, sep_val, true);
}

template <typename KdTraits>
void KdTree<KdTraits>::handle_extended_node(InternalPtr     nh,
                                            PointContainer &container,
                                            PointContainer &container_low)
{
	if (container_low.size() > 0)
	{
		nh->m_lower_low_val  = container_low.tbox().min_coord(nh->m_cut_dim);
		nh->m_lower_high_val = container_low.tbox().max_coord(nh->m_cut_dim);
	}
	else
	{
		nh->m_lower_low_val  = nh->m_cut_val;
		nh->m_lower_high_val = nh->m_cut_val;
	}
	if (container.size() > 0)
	{
		nh->m_upper_low_val  = container.tbox().min_coord(nh->m_cut_dim);
		nh->m_upper_high_val = container.tbox().max_coord(nh->m_cut_dim);
	}
	else
	{
		nh->m_upper_low_val  = nh->m_cut_val;
		nh->m_upper_high_val = nh->m_cut_val;
	}
}

} // namespace OMC