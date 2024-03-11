#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"
#include "OpenMeshCraft/Geometry/OrthogonalTree/OcTree.h"

#include "OpenMeshCraft/Mesh/IO/OBJReader.h"
#include "OpenMeshCraft/Mesh/IO/OBJWriter.h"
#include "OpenMeshCraft/Mesh/IO/STLReader.h"
#include "OpenMeshCraft/Mesh/TriSoup.h"

#include "OpenMeshCraft/Utils/Macros.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include "test_utils.h"

class test_OrthogonalTree : public testing::Test
{
protected:
	using index_t = OMC::index_t;
	using APAC    = OMC::APAC;

	class OrthogonalSplitPred
	{
	public:
		template <typename OrthTree, typename OrthNode>
		bool operator()(OMC_UNUSED const OrthTree &tree, const OrthNode &node)
		{
			return node.size() > 1000;
		}
	};

	class OrthogonalTraits
	{
	public:
		static constexpr size_t Dimension      = 3;
		static constexpr size_t MaxDepth       = 8;
		static constexpr bool   EnableVertices = false;

		using NT    = APAC::NT;
		using BboxT = APAC::BoundingBox3;

		struct VA
		{
			APAC::Point3 closest_p;
			double       dis;
		};
		using VertexAttrT = VA;

		using SplitPred   = OrthogonalSplitPred;
		using DoIntersect = APAC::DoIntersect;
		using CalcBbox    = APAC::CalcBoundingBox3;
	};

	using Tree = OMC::OcTree<OMC::OrthAutoDeduceTraits<OrthogonalTraits>>;

protected:
	Tree                    tree;
	static constexpr double dupl_thres = 2.5;

	Points    points;
	Triangles faces;

	void SetUp() override
	{
		TEST_GET_CONFIG(OrthogonalTree, SetUp);
		std::string filename = config.get<std::string>("filename");
		IOOptions   io_options;
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
		tree.construct(true, 1.5, dupl_thres, 64);
	}

	void TearDown() override { tree.clear(); }
};

TEST_F(test_OrthogonalTree, Construct)
{
	TEST_OUTPUT_DIRECTORY(OrthogonalTree, Construct);

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

		APAC::Point3 center           = tree.node_center(cur_node);
		APAC::Point3 half_side_length = tree.node_side_length(cur_node) * 0.5;

		//            xyz
		APAC::Point3 p000 = center - half_side_length;
		APAC::Point3 p111 = center + half_side_length;
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

		if (!cur_node.is_leaf())
		{
			// process each of its children
			for (index_t i = 0; i < Tree::Degree; ++i)
				nodes_to_traverse.push(cur_node.child(i));
		}
	}
	IOOptions io_options;
	io_options.vertex_has_point = true;
	write_mesh(outdir + "tree.obj", out_points, out_faces, io_options);

	std::fstream fout;
	fout.open(outdir + "dupl_tris.obj", std::ios::out);
	size_t fcnt = 0;
	for (index_t node_idx : tree.all_leaf_nodes())
	{
		if (tree.node(node_idx).dupl_degree() < dupl_thres)
			continue;

		for (const auto &box_ptr : tree.node(node_idx).boxes())
		{
			index_t     t_id = box_ptr->id();
			const auto &v0   = points[faces[t_id][0]];
			const auto &v1   = points[faces[t_id][1]];
			const auto &v2   = points[faces[t_id][2]];
			fout << "v " << v0.x() << " " << v0.y() << " " << v0.z() << "\n";
			fout << "v " << v1.x() << " " << v1.y() << " " << v1.z() << "\n";
			fout << "v " << v2.x() << " " << v2.y() << " " << v2.z() << "\n";
			fout << "f " << fcnt * 3 + 1 << " " << fcnt * 3 + 2 << " " << fcnt * 3 + 3
			     << "\n";
			fcnt++;
		}
	}
	fout.close();
}

#if 0
TEST_F(test_OrthogonalTree, Vertex)
{
	TEST_OUTPUT_DIRECTORY(OrthogonalTree, Construct);

	// find closest point for each vertex and save them as mesh to visualize.
	for (auto &v : tree.vertices())
	{
		v.attribute().dis = DBL_MAX;
	}

	std::queue<index_t> nodes_to_traverse;
	nodes_to_traverse.push(tree.root_node_idx());
	while (!nodes_to_traverse.empty())
	{
		index_t       cur_node_idx = nodes_to_traverse.front();
		Tree::NodeRef cur_node     = tree.node(cur_node_idx);
		nodes_to_traverse.pop();

		APAC::Point3 center = tree.node_center(cur_node);

		if (cur_node.is_leaf())
		{
			for (index_t i = 0; i < Tree::Degree; i++)
			{
				const APAC::Point3 &vp = tree.vertex(cur_node.vertex(i)).position();
				for (auto it : cur_node.boxes())
				{
					double dis = (points[it->id()] - vp).length();
					if (dis < tree.vertex(cur_node.vertex(i)).attribute().dis)
					{
						tree.vertex(cur_node.vertex(i)).attribute().dis = dis;
						tree.vertex(cur_node.vertex(i)).attribute().closest_p =
						  points[it->id()];
					}
				}
			}
		}

		if (!cur_node.is_leaf())
		{
			// process each of its children
			for (index_t i = 0; i < Tree::Degree; ++i)
				nodes_to_traverse.push(cur_node.child(i));
		}
	}

	IOOptions io_options;
	io_options.vertex_has_point = true;
	OBJWriter obj_writer;
	obj_writer.clear();
	for (auto v : tree.vertices())
	{
		if (v.attribute().dis < DBL_MAX)
		{
			index_t v0 = obj_writer.m_points.size();
			obj_writer.m_points.push_back(v.position());
			index_t v1 = obj_writer.m_points.size();
			obj_writer.m_points.push_back(v.attribute().closest_p);
			index_t v2 = obj_writer.m_points.size();
			obj_writer.m_points.push_back(
			  (v.position() + v.attribute().closest_p.as_vec()) * 0.5);
			obj_writer.m_triangles.push_back(Vec3id{v0, v1, v2});
		}
	}

	obj_writer.write(outdir + "vertices.obj", io_options, 10);
}

TEST_F(test_OrthogonalTree, LocatePoint)
{
	for (const auto &p : points)
	{
		index_t       node_idx = tree.locate(p);
		Tree::NodeRef node     = tree.node(node_idx);
		bool          found_p  = false;

		for (auto it : node.boxes())
		{
			if (points[it->id()] == p)
			{
				found_p = true;
				break;
			}
		}

		ASSERT_EQ(found_p, true);
	}
}
#endif