#pragma once

#include "OpenMeshCraft/Utils/IndexDef.h"
#include "OpenMeshCraft/Utils/InlinedVector.h"

namespace OMC {

/**
 * @brief Piecewise linear complex (PLC) for constraints recovery.
 */
template <typename Traits>
class PiecewiseLinearComplex
{
public: /* Traits **********************************************************/
  using Self = PiecewiseLinearComplex<Traits>;

  using NT      = typename Traits::NT;
  using EPoint3 = typename Traits::EPoint3;
  using GPoint3 = typename Traits::GPoint3;

  using AsGP = typename Traits::AsGP;
  using AsEP = typename Traits::AsEP;
  using ToEP = typename Traits::ToEP;

  using Orient3D         = typename Traits::Orient3D;
  using DotProductSign3D = typename Traits::DotProductSign3D;
  using InSphere         = typename Traits::InSphere;
  using CollinearPoints3 = typename Traits::CollinearPoints3;

public: /* Auxiliary data structures ****************************************/
  /// The type of the edge in the PLC, determining the recovery strategy on
  /// edges.
  enum class PLCEdgeType
  {
    // The type of the edge is undetermined.
    UNDETERMINED,
    // The "flat" edges will be ignored by segment recovery algorithm and will
    // not be further classyfied.
    FLAT_EDGE,
    // The "both acute vertex" edges will be immediatelly split by inserting the
    // middle point (each sub-edge becomes "one acute vertex")
    BOTH_ACUTE_VERTEX,
    // sub-edges of "one acute vertex" inherit the same edge type
    ONE_ACUTE_VERTEX,
    // sub-edges of "no acute vertex" inherit the same edge type
    NO_ACUTE_VERTEX,
    // The edge is marked to be deleted.
    TO_DELETE
  };

  class PLCEdge
  {
  public: /* Data *********************************************************/
    /// The type of the edge.
    PLCEdgeType type;
    /// The indices of the endpoints of the edge.
    /// (ep0 is always closer to the original ep0, same for ep1.
    /// So, the parameter `t` of ep0 is always less than `t` of ep1.)
    IPair       ep;
    /// the ancestor (not the parent) edge is the original edge before splitting
    /// (set to InvalidIndex if this edge is not a sub-edge)
    index_t     ancestor_id;
    /// the parent edge is the edge before splitting
    /// (set to InvalidIndex if this edge is not a sub-edge)
    index_t     parent_id;
    /// the child edge is the edge after splitting
    /// (set to InvalidIndex if this edge is not splitted)
    /// (every edge is splitted to two subsequent child edges)
    index_t     child_id;
    /// The original acute vertex is need to be remembered when splitting the
    /// edge with type `ONE_ACUTE_VERTEX`.
    index_t     acute_vid;

  public: /* Constructors *************************************************/
    PLCEdge() = default;
    PLCEdge(index_t e0, index_t e1);
    PLCEdge(PLCEdgeType _type, index_t e0, index_t e1, index_t _ancestor_id,
            index_t _parent_id, index_t _child_id, index_t _acute_vid);

    PLCEdge(const PLCEdge &rhs)            = default;
    PLCEdge(PLCEdge &&rhs)                 = default;
    PLCEdge &operator=(const PLCEdge &rhs) = default;
    PLCEdge &operator=(PLCEdge &&rhs)      = default;

  public: /* Interfaces ***************************************************/
    index_t       &ep0() { return ep.first; }
    const index_t &ep0() const { return ep.first; }
    index_t       &ep1() { return ep.second; }
    const index_t &ep1() const { return ep.second; }

    bool isDeleted() const { return type == PLCEdgeType::TO_DELETE; }
    bool isFlat() const { return type == PLCEdgeType::FLAT_EDGE; }

    bool hasAncestor() const { return is_valid_idx(ancestor_id); }
    bool hasParent() const { return is_valid_idx(parent_id); }
    bool isSplit() const { return is_valid_idx(child_id); }

    bool isConstraint() const
    {
      return !isDeleted() && !isFlat() && !isSplit();
    }

    bool hasEp(index_t vid) const { return ep0() == vid || ep1() == vid; }

    void makeUniqEp() { ep = unique_pair(ep.first, ep.second); }
    void swapEp() { std::swap(ep.first, ep.second); }

    bool operator<(const PLCEdge &rhs) const { return ep < rhs.ep; }
    bool operator==(const PLCEdge &rhs) const { return ep == rhs.ep; }

    index_t commonEp(const PLCEdge &rhs) const;
  };

  /// A range of sub-edges for an original edge after splitting.
  struct SubEdgeRange
  {
    index_t orig_eid; ///< the index of the original edge.
    index_t start; ///< the position in `sub_edges` where the sub-edges start.
    index_t size;  ///< the number of sub-edges stored in `sub_edges`.

    // clang-format off
    SubEdgeRange() = default;
    SubEdgeRange(index_t _eid, index_t _start, index_t _size)
      : orig_eid(_eid), start(_start), size(_size) {}
    SubEdgeRange(const SubEdgeRange &) = default;
    SubEdgeRange(SubEdgeRange &&)      = default;
    SubEdgeRange &operator=(const SubEdgeRange &rhs)
    { orig_eid = rhs.orig_eid; start = rhs.start; size = rhs.size; return *this; }
    // clang-format on
  };

  struct BoundingEdge
  {
    SubEdgeRange range;    // sub-edge range
    index_t      tid;      // adjacent triangle (not PLC face) index
    bool         reversed; // whether the edge is reversed

    // clang-format off
    bool operator==(index_t eid) const { return range.orig_eid == eid; }
    bool operator==(const BoundingEdge &rhs) const { return range.orig_eid == rhs.range.orig_eid; }
    bool operator!=(const BoundingEdge &rhs) const { return !(*this == rhs); }
    bool operator<(index_t eid) const { return range.orig_eid < eid; }
    bool operator<(const BoundingEdge &rhs) const { return range.orig_eid < rhs.range.orig_eid; }
    // clang-format on
  };

  class PLCFace
  {
  public: /* Data *********************************************************/
    /// Original triangles composing the face
    /// (triangle index, mutiple 3 with offset to point to `triangles`)
    InlinedVector2<index_t> triangles;

    /// Edges bounding the face (index to `plc_edges`)
    /// - The face may not be simply connected, i.e., there may be holes and
    /// multiple disjoint boundaries.
    /// - Duplicate bounding edges may exist during construction, these
    /// duplicate edges may be subsequent or not. We will remove the duplicate
    /// edges finally.
    InlinedVector4<BoundingEdge> bounding_edges;

    /// The unordered bounding vertices of the face (index to `vertices`)
    InlinedVector4<index_t> bounding_vertices;
    /// The unordered flat vertices of the face (index to `vertices`)
    InlinedVector4<index_t> flat_vertices;
  };

  // Inlined Vector Const Iterator
  using IVCIter = typename InlinedVector16<index_t>::const_iterator;

public: /* Constructor & Destructor ****************************************/
  PiecewiseLinearComplex() = delete;
  PiecewiseLinearComplex(const std::vector<GPoint3 *> &_vertices,
                         const std::vector<index_t>   &_edges,
                         const std::vector<index_t>   &_triangles);

public: /* Interfaces ******************************************************/
  /* Initialize PLC edges and faces ****************************************/
  void initPLCEdges();
  void initPLCFaces();

  /* Sub-parts of initializing PLC edges **********************************/

  void buildInitialPLCEdges();

  void classifyFlatEdges();

  void initVertIncEdge();

  void classifyVertEdge();

  /* Sub-parts of initializing PLC faces **********************************/

  void initSubEdges();

  void sortSubEdges();

  void assembleEdges2Faces();

  void mergeFacesArossFlatEdges();

  void removeDuplicateBoundingEdges();

  void extractBoundingVertices();

  /* Query on input constraints ***********************************************/

  index_t edgeVtx(index_t eid, index_t j) const { return edges[eid * 2 + j]; }

  index_t triVtx(index_t tid, index_t j) const
  {
    return triangles[tid * 3 + j];
  }

  /* Operations on PLC *******************************************************/

  /* Operations on PLC vertex */

  size_t         numVertices() const { return vertices.size(); }
  GPoint3       &pnt(index_t vid) { return *vertices[vid]; }
  const GPoint3 &pnt(index_t vid) const { return *vertices[vid]; }

  void newVtx(index_t new_vid);

  // clang-format off
  index_t& edgeWrtSteiner(index_t vid) { return edge_wrt_steiner[vid - input_nv]; }
  index_t edgeWrtSteiner(index_t vid) const { return edge_wrt_steiner[vid - input_nv]; }
  // clang-format on

  index_t steinerOfEdge(index_t eid) const;

  void removeIsolatedSteiner(index_t iso_vid);

  /* Operations on PLC edge */

  size_t         numEdges() const { return plc_edges.size(); }
  PLCEdge       &edge(index_t eid) { return plc_edges[eid]; }
  const PLCEdge &edge(index_t eid) const { return plc_edges[eid]; }

  index_t ancestorEdge(index_t eid) const;

  index_t parentEdge(index_t eid) const;

  index_t oppV2E(const PLCEdge &edge, index_t tid) const;

  void splitPLCEdge(index_t eid, index_t vid);

  index_t edgeExists(index_t ep0, index_t ep1) const;

  index_t collapsePLCEdge(index_t ep0, index_t ep1);

  void removeIsolatedEdges(index_t eid0, index_t eid1);

  /* Operations on PLC face */

  size_t         numFaces() const { return plc_faces.size(); }
  PLCFace       &face(index_t fid) { return plc_faces[fid]; }
  const PLCFace &face(index_t fid) const { return plc_faces[fid]; }

  /* Operations on PLC subedges */

  PLCEdge       &subEdge(index_t seid) { return edge(sub_edges[seid]); }
  const PLCEdge &subEdge(index_t seid) const { return edge(sub_edges[seid]); }

  /* Operations on PLC incidents (vert-edge, edge-tri, tri-bounding edge) */

  size_t  numEdgeIncTri(index_t eid) const { return edge_inc_tri[eid].size(); }
  index_t edgeIncTri(index_t eid, index_t j) { return edge_inc_tri[eid][j]; }

  std::pair<IVCIter, IVCIter> vertIncEdges(index_t vid) const;

  void updateVertIncEdge(index_t vid, index_t old_eid, index_t new_eid);

  const PLCEdge &boundingEdge(const PLCFace &f, index_t eid,
                              index_t *tid      = nullptr,
                              bool    *reversed = nullptr) const;

  template <typename IndexPairSet>
  void buildBoundingVtxAdjSet(const PLCFace &f, IndexPairSet &adj_vtx) const;

public: /* Data ************************************************************/
  /* Input data for PLC */

  /// The input vertices.
  const std::vector<GPoint3 *> &vertices;
  /// The input isolated edges (storing indices to vertices).
  const std::vector<index_t>   &edges;
  /// The input triangles (storing indices to vertices).
  const std::vector<index_t>   &triangles;

  const size_t input_nv; ///< number of input vertices
  const size_t input_ne; ///< number of input edges
  const size_t input_nt; ///< number of input triangles

  bool is_close_and_manifold; ///< whether the input is close and manifold

  /* Vertex data for PLC */

  /// Index to the edge with respect to its split Steiner point.
  /// (Steiner point index -> PLC edge index)
  std::vector<index_t> edge_wrt_steiner;

  /* Edge data for PLC */

  /// The built edges in the PLC.
  std::vector<PLCEdge> plc_edges;
  /// The initial number of PLC edges.
  size_t               init_npe;

  /* Face data for PLC */

  /// The built faces in the PLC.
  std::vector<PLCFace> plc_faces;

  /* Subedge data for PLC */

  /// Sub-edges range after splitting the original edges.
  std::vector<SubEdgeRange> sub_edge_range;
  /// All sub-edges are subsequently stored in this vector.
  std::vector<index_t>      sub_edges;

  /* Incidence data for PLC */

  /// Vertex --- Edge
  /// (We store incident edges for input vertices and split vertices separately)
  std::vector<InlinedVector16<index_t>> vertex_inc_edge_input;
  std::vector<InlinedVector2<index_t>>  vertex_inc_edge_steiner;

  /// Edge --- Triangle
  /// (We only store incident triangles for original input edges.)
  std::vector<InlinedVector4<index_t>> edge_inc_tri;
};

} // namespace OMC

#ifdef OMC_HAS_IMPL
  #include "PLC.inl"
#endif
