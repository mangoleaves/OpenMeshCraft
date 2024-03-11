#include "OpenMeshCraft/Geometry/AABB/AABBTree_Triangle_NearestSearch.h"
#include "OpenMeshCraft/Utils/StringUtils.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/AABB_traits.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_triangle_primitive.h"
#include "CGAL/Simple_cartesian.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

class test_AABBTreeClosestPoint : public testing::Test
{
protected:
	using APAC    = OMC::APAC;
	using OurTree = OMC::AABBTree_Triangle_NearestSearch;

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

protected:
	APAC::BoundingBox3          box;
	std::vector<APAC::Point3>   points;
	std::vector<OurTree::PrimT> triangles;
	std::vector<Triangle>       cgal_triangles;

protected:
	void SetUp() override
	{
		TEST_GET_CONFIG(AABBTreeClosestPoint, SetUp);

		std::ifstream fin;
		fin.open(config.get<std::string>("model_path"));
		ASSERT_EQ(fin.is_open(), true);

		std::string line;
		size_t      tri_id = 0;
		// read the mesh
		while (std::getline(fin, line))
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
				box += points.back();
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
				  OurTree::TriT(points[v0 - 1], points[v1 - 1], points[v2 - 1]),
				  tri_id);
				tri_id++;
			}
			break;
			}
		}
		for (auto &tri : triangles)
		{
			cgal_triangles.emplace_back(
			  Point(tri.v0().x(), tri.v0().y(), tri.v0().z()),
			  Point(tri.v1().x(), tri.v1().y(), tri.v1().z()),
			  Point(tri.v2().x(), tri.v2().y(), tri.v2().z()));
		}
	}

	bool double_equal(double dis, double cgal_dis)
	{
		if (std::abs(dis - cgal_dis) > 1e-10)
		{
			printf("ours %.15f, cgal %.15f\n", dis, cgal_dis);
		}
		return true;
	};
};

/**
 * @brief Compare closest point search by AABB with CGAL.
 */
TEST_F(test_AABBTreeClosestPoint, ClosestPoint)
{
	OurTree tree;
	// test three insert interfaces.
	tree.insert(triangles);
	tree.insert(triangles.begin(), triangles.end());
	tree.insert(std::move(triangles));
	// build and accelerate
	tree.build();
	tree.accelerate_nearest_search();

	// build cgal tree
	CGAL_Tree cgal_tree;
	cgal_tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	cgal_tree.build();
	cgal_tree.accelerate_distance_queries();

	double box_x = box.max_bound().x() - box.min_bound().x();
	double box_y = box.max_bound().y() - box.min_bound().y();
	double box_z = box.max_bound().z() - box.min_bound().z();

	double our_time = 0., cgal_time = 0.;

	// srand(time(nullptr));
	for (size_t i = 0; i < 100000; i++)
	{
		double x = (double)rand() / RAND_MAX * box_x;
		double y = (double)rand() / RAND_MAX * box_y;
		double z = (double)rand() / RAND_MAX * box_z;

		APAC::Point3 query(x, y, z);
		query += box.min_bound().as_vec();

		auto         start = OMC::Logger::elapse_reset();
		APAC::Point3 cp    = tree.closest_point(query);
		double       dis   = (cp - query).sqrnorm();
		our_time += OMC::Logger::elapsed(start).count();

		start = OMC::Logger::elapse_reset();
		Point  cgal_query(query.x(), query.y(), query.z());
		auto   cgal_result = cgal_tree.closest_point_and_primitive(cgal_query);
		double cgal_dis    = (cgal_query - cgal_result.first).squared_length();
		cgal_time += OMC::Logger::elapsed(start).count();

		ASSERT_TRUE(double_equal(dis, cgal_dis));
	}

	OMC::Logger::info(
	  std::format("our time {} s, CGAL time {} s\n", our_time, cgal_time));
}