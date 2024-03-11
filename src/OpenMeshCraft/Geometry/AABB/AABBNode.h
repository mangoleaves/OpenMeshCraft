#pragma once

namespace OMC {

/**
 * @brief Define node of AABB tree.
 * @tparam PrimT The geometry primitive type stored in leaf nodes.
 * @tparam BboxT The bounding box type used in internal nodes.
 */
template <typename PrimT, typename BboxT>
class AABBNode
{
public:
	using NodeT = AABBNode<PrimT, BboxT>;

	/// pointer to primitive
	using PrimPtr   = PrimT *;
	using PrimCPtr  = const PrimT *;
	/// pointer to AABBNode
	using NodePtr   = NodeT *;
	using NodeCPtr  = const NodeT *;
	/// a general pointer that can be reinterpreted as first two pointers
	using ChildPtr  = void *;
	using ChildCPtr = const void *;

public:
	AABBNode()
	  : m_bbox()
	  , m_p_left_child(nullptr)
	  , m_p_right_child(nullptr)
	{
	}

	~AABBNode() {}

	inline BboxT       &bbox() { return m_bbox; }
	inline const BboxT &bbox() const { return m_bbox; }

	inline ChildPtr        &left_ptr() { return m_p_left_child; }
	inline const ChildCPtr &left_ptr() const { return m_p_left_child; }
	inline ChildPtr        &right_ptr() { return m_p_right_child; }
	inline const ChildCPtr &right_ptr() const { return m_p_right_child; }

	inline NodePtr left_child()
	{
		return reinterpret_cast<NodePtr>(m_p_left_child);
	}
	inline NodePtr right_child()
	{
		return reinterpret_cast<NodePtr>(m_p_right_child);
	}

	inline PrimPtr left_data()
	{
		return reinterpret_cast<PrimCPtr>(m_p_left_child);
	}
	inline PrimPtr right_data()
	{
		return reinterpret_cast<PrimCPtr>(m_p_right_child);
	}

	inline NodeCPtr left_child() const
	{
		return reinterpret_cast<NodeCPtr>(m_p_left_child);
	}
	inline NodeCPtr right_child() const
	{
		return reinterpret_cast<NodeCPtr>(m_p_right_child);
	}

	inline PrimCPtr left_data() const
	{
		return reinterpret_cast<PrimCPtr>(m_p_left_child);
	}
	inline PrimCPtr right_data() const
	{
		return reinterpret_cast<PrimCPtr>(m_p_right_child);
	}

private:
	BboxT    m_bbox;
	ChildPtr m_p_left_child, m_p_right_child;
};

} // namespace OMC
