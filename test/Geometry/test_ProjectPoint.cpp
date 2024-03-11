#include "OpenMeshCraft/Geometry/ApproxPredicatesApproxConstructions.h"

#include "AreVecsEqual.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Simple_cartesian.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

#include <type_traits>

class test_ProjectPoint : public testing::Test
{
protected:
	using Our_K    = OMC::APAC;
	using Our_Pnt3 = Our_K::Point3;
	using Our_Tri3 = Our_K::Triangle3;
	class Our_Htri3 : public Our_Tri3,
	                  public Our_K::FastProjectPoint3::AuxTriangle
	{
	public:
		using NT           = typename Our_K::NT;
		using VecT         = typename Our_K::Vec3;
		using PointT       = typename Our_K::Point3;
		using BaseTriangle = typename Our_K::Triangle3;
		using AuxTriangle  = typename Our_K::FastProjectPoint3::AuxTriangle;

	public:
		Our_Htri3(const PointT &v0, const PointT &v1, const PointT &v2)
		  : BaseTriangle(v0, v1, v2)
		  , AuxTriangle(v0, v1, v2)
		{
		}
	};
	using Our_ProjectPoint     = Our_K::ProjectPoint3;
	using Our_FastProjectPoint = Our_K::FastProjectPoint3;

	using CGAL_K    = CGAL::Simple_cartesian<double>;
	using CGAL_Pnt3 = CGAL::Point_3<CGAL_K>;
	using CGAL_Tri3 = CGAL::Triangle_3<CGAL_K>;

	bool pnts_equal(const Our_Pnt3 &t, const CGAL_Pnt3 &s)
	{
		for (size_t i = 0; i < 3; i++)
		{
			if (std::abs(t[i] - s[(int)i]) > 1e-12)
			{
				std::cout << std::format("{}-th elements are not equal: {} in != {}.",
				                         i, t[i], s[(int)i]);
				return true;
			}
		}
		return true;
	};
};

TEST_F(test_ProjectPoint, Projection)
{
	srand((unsigned int)time(nullptr));

	CGAL_Pnt3 p0{1.11, 2.22, 3.33}, p1{3.14159, 2.7, 1.414}, p2{5., 6., 7.};
	Our_Pnt3  v0{1.11, 2.22, 3.33}, v1{3.14159, 2.7, 1.414}, v2{5., 6., 7.};

	CGAL_Tri3 cgal_tri(p0, p1, p2);
	Our_Tri3  our_tri(v0, v1, v2);
	Our_Htri3 our_htri(v0, v1, v2);

	CGAL_Pnt3 query_p;
	Our_Pnt3  query_v;

	auto   our_projct     = Our_ProjectPoint();
	auto   our_fastprojct = Our_FastProjectPoint();
	CGAL_K kernel;
	auto   project = kernel.construct_projected_point_3_object();

	for (size_t i = 0; i < 1000; i++)
	{
		query_v.x() = ((double)rand() / RAND_MAX) * 20. - 5.0;
		query_v.y() = (double)rand() / RAND_MAX * 20. - 5.0;
		query_v.z() = (double)rand() / RAND_MAX * 20. - 5.0;

		query_p = CGAL_Pnt3(query_v.x(), query_v.y(), query_v.z());

		CGAL_Pnt3 project_p  = project(cgal_tri, query_p);
		Our_Pnt3  project_v  = our_projct(our_tri, query_v);
		Our_Pnt3  project_hv = our_fastprojct(our_htri, query_v);

		ASSERT_TRUE(pnts_equal(project_v, project_p));
		ASSERT_TRUE(pnts_equal(project_hv, project_p));
	}
}