#pragma once

namespace OMC {

/**
 * @brief Common node type for KdTree.
 */
template <typename KdTraits>
class KdNode
{
public:
	KdNode() = delete;

	KdNode(bool leaf)
	  : m_is_leaf(leaf)
	{
	}

public:
	bool m_is_leaf;
};

/**
 * @brief Internal node type for KdTree.
 */
template <typename KdTraits>
class KdInternalNode : public KdNode<KdTraits>
{
public:
	using NT      = typename KdTraits::NT;
	using Node    = KdNode<KdTraits>;
	using NodePtr = Node *;

public:
	KdInternalNode()
	  : Node(false)
	{
	}

	inline void set_separator(size_t cut_dim, NT cut_val)
	{
		m_cut_dim = cut_dim;
		m_cut_val = cut_val;
	}

public:
	size_t  m_cut_dim;
	NT      m_cut_val;
	NT      m_lower_low_val, m_lower_high_val, m_upper_low_val, m_upper_high_val;
	NodePtr m_lower_ch, m_upper_ch;
};

/**
 * @brief Leaf node type of KdTree.
 */
template <typename KdTraits>
class KdLeafNode : public KdNode<KdTraits>
{
public:
	using NT           = typename KdTraits::NT;
	using Node         = KdNode<KdTraits>;
	using KdPointsIter = typename KdTraits::KdPointsIter;

public:
	KdLeafNode()
	  : Node(true)
	{
	}

	KdLeafNode(size_t n)
	  : Node(true)
	  , n(n)
	{
	}

public:
	size_t       n;
	KdPointsIter data;
};

} // namespace OMC