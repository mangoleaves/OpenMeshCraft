#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/BinaryTree/BinaryTree.h"
#include "OpenMeshCraft/Geometry/BinaryTree/BinarySplitManner.h"

#include "OpenMeshCraft/Utils/Macros.h"

#include "test_utils.h"

class test_BinaryTree : public testing::Test
{
protected:
	using index_t = OMC::index_t;
	using APAC    = OMC::APAC;

	class BinarySplitPred
	{
	public:
		template <typename BinaryTree, typename BinaryNode>
		bool operator()(OMC_UNUSED const BinaryTree &tree, const BinaryNode &node)
		{
			return node.size() > 50;
		}
	};

	class BinaryShapeRefinePred
	{
	public:
		template <typename BinaryTree, typename BinaryNode>
		bool operator()(const BinaryTree &tree, const BinaryNode &node,
		                std::array<bool, 3> &partitionable)
		{
			auto diag_length =
			  (tree.box().max_bound() - tree.box().min_bound()).length();
			auto node_length = node.box().max_bound() - node.box().min_bound();
			partitionable[0] = node_length[0] > 0.05 * diag_length;
			partitionable[1] = node_length[1] > 0.05 * diag_length;
			partitionable[2] = node_length[2] > 0.05 * diag_length;
			return partitionable[0] || partitionable[1] || partitionable[2];
		}
	};

	class BinaryTraits
	{
	public:
		static constexpr size_t Dimension = 3;
		static constexpr size_t MaxDepth  = 32;

		using NT    = APAC::NT;
		using BboxT = APAC::BoundingBox3;

		using SplitPred       = BinarySplitPred;
		using SplitManner     = OMC::BinarySplitManner;
		using ShapeRefinePred = BinaryShapeRefinePred;
		using DoIntersect     = APAC::DoIntersect;
		using CalcBbox        = APAC::CalcBoundingBox3;
	};

	using Tree = OMC::BinaryTree<OMC::BinaryAutoDeduceTraits<BinaryTraits>>;

protected:
	Tree tree;

	Points    points;
	Triangles faces;

	void SetUp() override
	{
		TEST_GET_CONFIG(BinaryTree, SetUp);
		std::string filename = config.get<std::string>("filename");

		IOOptions io_options;
		io_options.vertex_has_point = true;

		read_mesh(filename, points, faces, io_options);

		std::vector<index_t>         indices;
		std::vector<APAC::Triangle3> triangles;

		for (const auto &f : faces)
		{
			triangles.push_back(
			  APAC::Triangle3(points[f[0]], points[f[1]], points[f[2]]));
			indices.push_back(indices.size());
		}

		tree.insert_primitives(triangles, indices);
		tree.construct(true, 1.01);
		// tree.shape_refine();
	}

	void TearDown() override { tree.clear(); }
};

TEST_F(test_BinaryTree, Construct)
{
	TEST_OUTPUT_DIRECTORY(BinaryTree, Construct);

	// visualize tree and save
	Points    out_points;
	Triangles out_faces;

	std::queue<index_t> nodes_to_traverse;
	nodes_to_traverse.push(tree.root_node_idx());

	out_points.clear();
	out_faces.clear();

	while (!nodes_to_traverse.empty())
	{
		index_t        cur_node_idx = nodes_to_traverse.front();
		Tree::NodeCRef cur_node     = tree.node(cur_node_idx);
		nodes_to_traverse.pop();

		if (!cur_node.is_leaf())
		{
			// process each of its children
			for (index_t i = 0; i < cur_node.children_size(); ++i)
				nodes_to_traverse.push(cur_node.child(i));
			continue;
		}

		//            xyz
		APAC::Point3 p000 = cur_node.box().min_bound();
		APAC::Point3 p111 = cur_node.box().max_bound();
		APAC::Point3 p001 = p000, p010 = p000, p011 = p111, p100 = p000,
		             p101 = p111, p110 = p111;
		p001.z() = p111.z();
		p010.y() = p111.y();
		p011.x() = p000.x();
		p100.x() = p111.x();
		p101.y() = p000.y();
		p110.z() = p000.z();

		index_t v000 = out_points.size();
		out_points.push_back(p000);
		index_t v001 = out_points.size();
		out_points.push_back(p001);
		index_t v010 = out_points.size();
		out_points.push_back(p010);
		index_t v011 = out_points.size();
		out_points.push_back(p011);
		index_t v100 = out_points.size();
		out_points.push_back(p100);
		index_t v101 = out_points.size();
		out_points.push_back(p101);
		index_t v110 = out_points.size();
		out_points.push_back(p110);
		index_t v111 = out_points.size();
		out_points.push_back(p111);

		out_faces.emplace_back(v000, v100, out_points.size());
		out_points.push_back((p000 + p100.as_vec()) * 0.5);
		out_faces.emplace_back(v000, v010, out_points.size());
		out_points.push_back((p000 + p010.as_vec()) * 0.5);
		out_faces.emplace_back(v000, v001, out_points.size());
		out_points.push_back((p000 + p001.as_vec()) * 0.5);
		out_faces.emplace_back(v100, v110, out_points.size());
		out_points.push_back((p100 + p110.as_vec()) * 0.5);
		out_faces.emplace_back(v100, v101, out_points.size());
		out_points.push_back((p100 + p101.as_vec()) * 0.5);
		out_faces.emplace_back(v010, v110, out_points.size());
		out_points.push_back((p010 + p110.as_vec()) * 0.5);
		out_faces.emplace_back(v010, v011, out_points.size());
		out_points.push_back((p010 + p011.as_vec()) * 0.5);
		out_faces.emplace_back(v110, v111, out_points.size());
		out_points.push_back((p110 + p111.as_vec()) * 0.5);
		out_faces.emplace_back(v001, v101, out_points.size());
		out_points.push_back((p001 + p101.as_vec()) * 0.5);
		out_faces.emplace_back(v001, v011, out_points.size());
		out_points.push_back((p001 + p011.as_vec()) * 0.5);
		out_faces.emplace_back(v101, v111, out_points.size());
		out_points.push_back((p101 + p111.as_vec()) * 0.5);
		out_faces.emplace_back(v011, v111, out_points.size());
		out_points.push_back((p011 + p111.as_vec()) * 0.5);
	}
	IOOptions io_options;
	io_options.vertex_has_point = true;
	write_mesh(outdir + "tree.obj", out_points, out_faces, io_options);
}