#include "OpenMeshCraft/Geometry/AABB/AABBTree_Triangle_Intersection.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "CGAL/Simple_cartesian.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

class test_AABBTreeIntersection : public testing::Test
{
protected:
	// Our
	using APAC   = OMC::APAC;
	using MyTree = OMC::AABBTree_Triangle_Intersection<APAC>;
	using IndexedTriangle =
	  OMC::PrimitiveWithAttribute<APAC::Triangle3, size_t>;

	// CGAL
	using K         = typename CGAL::Simple_cartesian<double>;
	using FT        = typename K::FT;
	using Ray       = typename K::Ray_3;
	using Line      = typename K::Line_3;
	using Point     = typename K::Point_3;
	using Triangle  = typename K::Triangle_3;
	using Iterator  = typename std::vector<Triangle>::iterator;
	using Primitive = typename CGAL::AABB_triangle_primitive<K, Iterator>;
	using AABB_triangle_traits = typename CGAL::AABB_traits<K, Primitive>;
	using CGAL_Tree            = typename CGAL::AABB_tree<AABB_triangle_traits>;
	using Primitive_id         = typename CGAL_Tree::Primitive_id;

protected:
	std::vector<APAC::Point3>    points;
	std::vector<IndexedTriangle> triangles;
	std::vector<Triangle>        cgal_triangles;

protected:
	void SetUp() override
	{
		TEST_GET_CONFIG(AABBTreeIntersection, SetUp);

		std::ifstream ifs;
		ifs.open(config.get<std::string>("model_path"));
		EXPECT_EQ(ifs.is_open(), true);
		std::string line;
		size_t      tri_id = 0;
		// read the mesh
		while (std::getline(ifs, line))
		{
			switch (line[0])
			{
			case 'v':
			{
				auto   sub_strs = OMC::split_string(line, ' ');
				double x        = std::stod(sub_strs[1]);
				double y        = std::stod(sub_strs[2]);
				double z        = std::stod(sub_strs[3]);
				points.emplace_back(x, y, z);
			}
			break;
			case 'f':
			{
				// get face vertices
				auto sub_strs = OMC::split_string(line, ' ');
				int  v0       = std::stoi(sub_strs[1]);
				int  v1       = std::stoi(sub_strs[2]);
				int  v2       = std::stoi(sub_strs[3]);
				triangles.emplace_back(
				  APAC::Triangle3(points[v0 - 1], points[v1 - 1], points[v2 - 1]),
				  tri_id);
				tri_id++;
			}
			break;
			}
		}
		ifs.close();

		// construct CGAL triangles
		for (auto &tri : triangles)
		{
			cgal_triangles.emplace_back(
			  Point(tri.v0().x(), tri.v0().y(), tri.v0().z()),
			  Point(tri.v1().x(), tri.v1().y(), tri.v1().z()),
			  Point(tri.v2().x(), tri.v2().y(), tri.v2().z()));
		}
	}
};

TEST_F(test_AABBTreeIntersection, BoxIntersect)
{
	TEST_OUTPUT_DIRECTORY(AABBTreeIntersection, BoxIntersect);

	// construct our tree

	MyTree tree;
	tree.insert(triangles);
	tree.build();

	// test our intersection

	auto start = OMC::Logger::elapse_reset();

	std::vector<size_t> results;
	APAC::Point3        center(0.0, 0.0, 0.0);
	for (int j = 0; j < 3; ++j)
	{
		center.x() += triangles[0][j].x();
		center.y() += triangles[0][j].y();
		center.z() += triangles[0][j].z();
	}
	center.x() /= 3.0;
	center.y() /= 3.0;
	center.z() /= 3.0;

	APAC::Point3   aux_point1(points[0][0] + 0.01, points[0][1], points[0][2]);
	APAC::Point3   aux_point2(center[0] - 0.01, center[1], center[2]);
	APAC::Segment3 seg(aux_point1, aux_point2);
	tree.all_intersections(seg, results);

	std::cout << std::format("OpenMeshCraft AABBTree time: {} s\n",
	                         OMC::Logger::elapsed(start).count());

	// construct CGAL tree

	CGAL_Tree cgal_tree;
	cgal_tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	cgal_tree.build();

	// test CGAL intersection

	start = OMC::Logger::elapse_reset();

	std::list<Primitive_id> primitives;
	cgal_tree.all_intersected_primitives(cgal_triangles[1],
	                                     std::back_inserter(primitives));

	std::cout << std::format("CGAL AABBTree time: {} s\n",
	                         OMC::Logger::elapsed(start).count());

	// output intersection results for two methods

	std::ofstream ofs;

	make_file_writable(outdir + "cgal_intersect.txt");
	ofs.open(outdir + "OpenMeshCraft_intersect.txt");

	ofs << aux_point1[0] << " " << aux_point1[1] << " " << aux_point1[2]
	    << std::endl;

	ofs << aux_point2[0] << " " << aux_point2[1] << " " << aux_point2[2]
	    << std::endl;

	ofs << results.size() << std::endl;
	for (size_t i = 0; i < results.size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			ofs << triangles[results[i]][j].x() << " " << triangles[results[i]][j].y()
			    << " " << triangles[results[i]][j].z() << std::endl;
		}
	}
	ofs.close();

	make_file_writable(outdir + "cgal_intersect.txt");
	ofs.open(outdir + "cgal_intersect.txt");
	ofs << primitives.size() << std::endl;
	for (std::list<Primitive_id>::iterator it = primitives.begin();
	     it != primitives.end(); ++it)
	{
		const Primitive_id &primitive_id = *it;
		for (int i = 0; i < 3; ++i)
			ofs << primitive_id->vertex(i) << std::endl;
	}
	ofs.close();
}