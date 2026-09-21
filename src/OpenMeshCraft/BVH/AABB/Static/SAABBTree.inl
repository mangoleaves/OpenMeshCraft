#pragma once

#include "SAABBTree.h"

#include "tbb/tbb.h"

#include "OpenMeshCraft/Utils/ExecutionCompat.h"

namespace OMC {

template <typename Traits>
SAABBTree<Traits>::SAABBTree(Prims &&primitives)
{
  insert(std::move(primitives));
  build();
}

template <typename Traits>
template <typename PrimsIterT>
SAABBTree<Traits>::SAABBTree(PrimsIterT first, PrimsIterT beyond)
{
  insert(first, beyond);
  build();
}

template <typename Traits>
SAABBTree<Traits>::~SAABBTree()
{
  clear();
}

template <typename Traits>
void SAABBTree<Traits>::insert(Prims &&primitives)
{
  clear();
  m_primitives = std::move(primitives);
}

template <typename Traits>
template <typename PrimsIterT>
void SAABBTree<Traits>::insert(PrimsIterT first, PrimsIterT beyond)
{
  clear();
  m_primitives.insert(m_primitives.end(), first, beyond);
}

template <typename Traits>
void SAABBTree<Traits>::build(bool parallel)
{
  // calculate bounding box for each primitive
  update_boxes(parallel);

  // clear nodes and avoid wasting memory
  m_nodes.clear();
  m_nodes.shrink_to_fit();

  // If the number of primitives is less than 2, don't build the tree.
  if (m_primitives.size() <= 1)
    return;

  // we gonna construct a certain number of nodes
  m_nodes.resize(m_primitives.size() - 1);

  // sequentially split nodes until there is no node to split
  // or we get enough nodes to split them parallelly.

  const size_t NUM_NODES_FOR_PARALLEL =
    parallel ? (size_t)tbb::this_task_arena::max_concurrency()
             : std::numeric_limits<size_t>::max();

  Indices prim_indices(m_primitives.size());
  std::iota(prim_indices.begin(), prim_indices.end(), 0);

  struct NodeToSplit
  {
    index_t     node_id;
    IndicesIter first, beyond;
    // clang-format off
    NodeToSplit(index_t _node_id, IndicesIter _first, IndicesIter _beyond)
      : node_id(_node_id) , first(_first) , beyond(_beyond) {}
    // clang-format on
  };

  std::deque<NodeToSplit> nodes_to_split;
  nodes_to_split.emplace_back(/*root_node_id*/ 0, prim_indices.begin(),
                              prim_indices.end());

  while (!nodes_to_split.empty() &&
         nodes_to_split.size() < NUM_NODES_FOR_PARALLEL)
  {
    const NodeToSplit &cur_node_to_split = nodes_to_split.front();

    index_t cur_node_id = cur_node_to_split.node_id;
    NodeRef cur_node    = node(cur_node_id);

    IndicesIter first  = cur_node_to_split.first,
                beyond = cur_node_to_split.beyond;
    IndicesIter middle =
      expand(cur_node_id, first, beyond, parallel, /*recursive*/ false);

    size_t range = beyond - first;
    switch (range)
    {
    case 2:
      break;
    case 3:
      nodes_to_split.emplace_back(cur_node.right_child(), middle, beyond);
      break;
    default:
      nodes_to_split.emplace_back(cur_node.left_child(), first, middle);
      nodes_to_split.emplace_back(cur_node.right_child(), middle, beyond);
      break;
    }
    nodes_to_split.pop_front();
  }

  // if there remain enough nodes, split them parallelly.
  if (!nodes_to_split.empty())
  {
    tbb::parallel_for_each(nodes_to_split.begin(), nodes_to_split.end(),
                           [this](const NodeToSplit &root_node)
                           {
                             expand(root_node.node_id, root_node.first,
                                    root_node.beyond,
                                    /*parallel*/ false, /*recursive*/ true);
                           });
  }
}

template <typename Traits>
void SAABBTree<Traits>::update(bool parallel)
{
  update_boxes(parallel);

  // there is no node to update
  if (size() <= 1)
    return;

  size_t num_threads = (size_t)tbb::this_task_arena::max_concurrency();

  // if the number of threads is small or the number of primitives is small,
  // we can update the tree sequentially.
  if (!parallel || num_threads <= 1 || size() <= num_threads * num_threads ||
      size() < 64)
  {
    update_recursive(/*root_node_id*/ 0, m_primitives.size());
    return;
  }

  // update the tree parallelly
  size_t parallel_tasks = 1;
  while (parallel_tasks < num_threads)
    parallel_tasks *= 2;

  if (parallel_tasks >= size() / 4)
  {
    update_recursive(/*root_node_id*/ 0, m_primitives.size());
    return;
  }

  struct NodeToUpdate
  {
    index_t node_id;
    size_t  nb_primitives;
    // clang-format off
    NodeToUpdate(index_t _node_id, size_t _nb_primitives)
      : node_id(_node_id) , nb_primitives(_nb_primitives) {}
    // clang-format on
  };

  std::deque<NodeToUpdate> nodes_to_update;
  nodes_to_update.emplace_back(/*root_node_id*/ 0, m_primitives.size());

  // width-first-search nodes, prepare for parallel updating
  while (nodes_to_update.size() < parallel_tasks)
  {
    const NodeToUpdate &cur_node_to_split = nodes_to_update.front();

    index_t node_id       = cur_node_to_split.node_id;
    NodeRef cur_node      = node(node_id);
    size_t  nb_primitives = cur_node_to_split.nb_primitives;

    index_t left_node_id  = cur_node.left_child(),
            right_node_id = cur_node.right_child();

    size_t left_primitives  = nb_primitives / 2;
    size_t right_primitives = nb_primitives - left_primitives;

    nodes_to_update.emplace_back(left_node_id, left_primitives);
    nodes_to_update.emplace_back(right_node_id, right_primitives);

    nodes_to_update.pop_front();
  }

  // parallel updating
  tbb::parallel_for_each(nodes_to_update.begin(), nodes_to_update.end(),
                         [this](const NodeToUpdate &node_to_update)
                         {
                           update_recursive(node_to_update.node_id,
                                            node_to_update.nb_primitives);
                         });

  // width-first-search, find remain nodes to update sequentially
  std::deque<index_t> remain_nodes;
  remain_nodes.push_back(/*root_node_id*/ 0);
  size_t front = 0;

  while (remain_nodes.size() < parallel_tasks)
  {
    index_t front_node_id = remain_nodes[front];
    remain_nodes.push_back(node(front_node_id).left_child());
    remain_nodes.push_back(node(front_node_id).right_child());
    front += 1;
  }
  // update remain nodes sequentially
  for (int i = (int)remain_nodes.size() - 1; i >= 0; i--)
  {
    node(remain_nodes[i]).bbox() =
      node(node(remain_nodes[i]).left_child()).bbox() +
      node(node(remain_nodes[i]).right_child()).bbox();
  }
}

template <typename Traits>
void SAABBTree<Traits>::clear()
{
  m_nodes.clear();
  m_primitives.clear();
}

template <typename Traits>
template <typename TraversalTrait>
void SAABBTree<Traits>::traverse(TraversalTrait &traits) const
{
  switch (size())
  {
  case 0:
    break;
  case 1:
    if (traits.do_inter(m_prim_boxes[0]))
      traits.intersection(m_primitives[0]);
    break;
  default:
    traverse_node(/*root_node_id*/ 0, traits, m_primitives.size());
  }
}

template <typename Traits>
void SAABBTree<Traits>::update_boxes(bool parallel)
{
  m_prim_boxes.resize(m_primitives.size());
  if (parallel)
    std::transform(OMC_IF_EXEC(std::execution::par_unseq) m_primitives.begin(),
                   m_primitives.end(), m_prim_boxes.begin(), CalcBbox());
  else
    std::transform(m_primitives.begin(), m_primitives.end(),
                   m_prim_boxes.begin(), CalcBbox());
}

/**
 * @brief Expand the node to two sub-nodes.
 * @param node_id The node to be expanded.
 * @param first The first iterator of the primitives.
 * @param beyond The beyond iterator of the primitives.
 * @param parallel Whether to expand the node parallelly.
 * @param recursive Whether to expand the node recursively.
 * @return The middle iterator of the separated primitives.
 */
template <typename Traits>
auto SAABBTree<Traits>::expand(index_t node_id, IndicesIter first,
                               IndicesIter beyond, bool parallel,
                               bool recursive) -> IndicesIter
{
  NodeRef cur_node = node(node_id);

  // Calculate box for primitves(first, beyond).
  class ReduceOp
  {
    // clang-format off
  public:
    const SAABBTree<Traits> *tree;
    ReduceOp(const SAABBTree<Traits> *_tree) : tree(_tree) {}
  public:
    BboxT operator()(const BboxT &lhs_box, index_t rhs_prim_id) const { return lhs_box + tree->bbox(rhs_prim_id); }
    BboxT operator()(index_t lhs_prim_id, const BboxT &rhs_box) const { return tree->bbox(lhs_prim_id) + rhs_box; }
    BboxT operator()(const BboxT &lhs_box, const BboxT &rhs_box) const { return lhs_box + rhs_box; }
    BboxT operator()(index_t lhs_prim_id, index_t rhs_prim_id) const { return tree->bbox(lhs_prim_id) + tree->bbox(rhs_prim_id); }
    // clang-format on
  };
  BboxT init_box = bbox(*first);

#if defined(__cpp_lib_parallel_algorithm) && \
  __cpp_lib_parallel_algorithm >= 201603L
  if (parallel)
    cur_node.bbox() =
      std::reduce(OMC_IF_EXEC(std::execution::par_unseq) first + 1, beyond,
                  init_box, ReduceOp(this));
  else
    cur_node.bbox() = std::reduce(first + 1, beyond, init_box, ReduceOp(this));
#else
  cur_node.bbox() =
    std::accumulate(first + 1, beyond, init_box, ReduceOp(this));
#endif

  // Split primitives to two groups.
  IndicesIter middle;
  split_primitives(first, beyond, middle, cur_node.bbox(), parallel);

  size_t range = beyond - first;
  switch (range)
  {
  case 2:
    cur_node.left_child()  = *first;
    cur_node.right_child() = *(first + 1);
    // stop expanding
    return beyond;
  case 3:
    cur_node.left_child()  = *first;
    cur_node.right_child() = node_id + 1;
    if (recursive)
      expand(cur_node.right_child(), first + 1, beyond, parallel, recursive);
    return first + 1;
  default:
    size_t new_range       = middle - first;
    cur_node.left_child()  = node_id + 1;
    cur_node.right_child() = node_id + new_range;
    if (recursive)
    {
      expand(cur_node.left_child(), first, middle, parallel, recursive);
      expand(cur_node.right_child(), middle, beyond, parallel, recursive);
    }
    return middle;
  }
}

/**
 * @brief split a range of primitives based on a specified splitting predicate
 * and the longest axis of the bounding box. The primitives are always split at
 * the middle.
 * @param [in] first The first iterator of the primitives.
 * @param [in] beyond The beyond iterator of the primitives.
 * @param [out] middle The middle iterator of the primitives.
 * @param [in] box The bounding box of the primitives.
 * @return The longest axis (split axis) of the bounding box.
 */
template <typename Traits>
size_t SAABBTree<Traits>::split_primitives(IndicesIter  first,
                                           IndicesIter  beyond,
                                           IndicesIter &middle,
                                           const BboxT &box, bool parallel)
{
  size_t        split_axis = box.longest_axis();
  PrimSplitPred pred(split_axis);

  auto wrap_pred = [this, &pred](index_t lhs, index_t rhs)
  { return pred(primitive(lhs), primitive(rhs)); };

  middle = first + (beyond - first) / 2;
  if (parallel)
    std::nth_element(OMC_IF_EXEC(std::execution::par_unseq) first, middle,
                     beyond, wrap_pred);
  else
    std::nth_element(first, middle, beyond, wrap_pred);
  return split_axis;
}

template <typename Traits>
void SAABBTree<Traits>::update_recursive(index_t node_id, size_t nb_primitives)
{
  NodeRef nd = node(node_id);
  switch (nb_primitives)
  {
  case 2:
    nd.bbox() = bbox(nd.left_child()) + bbox(nd.right_child());
    break;
  case 3:
    update_recursive(nd.right_child(), 2);
    nd.bbox() = bbox(nd.left_child()) + node(nd.right_child()).bbox();
    break;
  default:
    update_recursive(nd.left_child(), nb_primitives / 2);
    update_recursive(nd.right_child(), nb_primitives - nb_primitives / 2);
    nd.bbox() = node(nd.left_child()).bbox() + node(nd.right_child()).bbox();
    break;
  }
}

template <typename Traits>
template <typename TraversalTrait>
bool SAABBTree<Traits>::traverse_node(index_t node_id, TraversalTrait &trait,
                                      size_t nb_primitives) const
{
  NodeCRef nd      = node(node_id);
  bool     go_next = true;
  switch (nb_primitives)
  {
  case 2:
    go_next = trait.intersection(primitive(nd.left_child()));
    if (go_next)
      go_next = trait.intersection(primitive(nd.right_child()));
    return go_next;
  case 3:
    go_next = trait.intersection(primitive(nd.left_child()));
    if (go_next && trait.do_inter(node(nd.right_child()).bbox()))
      go_next = traverse_node(nd.right_child(), trait, 2);
    return go_next;
  default:
    if (trait.do_inter(node(nd.left_child()).bbox()))
    {
      go_next = traverse_node(nd.left_child(), trait, nb_primitives / 2);
      if (go_next && trait.do_inter(node(nd.right_child()).bbox()))
        go_next = traverse_node(nd.right_child(), trait,
                                nb_primitives - nb_primitives / 2);
      return go_next;
    }
    else if (trait.do_inter(node(nd.right_child()).bbox()))
    {
      return traverse_node(nd.right_child(), trait,
                           nb_primitives - nb_primitives / 2);
    }
    return true;
  }
}

} // namespace OMC
