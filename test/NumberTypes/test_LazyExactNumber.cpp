#include "OpenMeshCraft/NumberTypes/LazyNumber.h"

#include "OpenMeshCraft/Utils/DisableWarnings.h"

#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Lazy_exact_nt.h"

#include "OpenMeshCraft/Utils/EnableWarnings.h"

#include "test_utils.h"

#include <type_traits>

using namespace OMC;

using Lazy_exact_nt = CGAL::Lazy_exact_nt<OMC::BoostRational>;

#define LazyNumber LazyNumber<std::true_type, std::true_type, std::true_type>

bool interval_number_equal(const IntervalNumber<std::true_type> &it,
                           const CGAL::Interval_nt<false>       &_it)
{
	if (!(it.inf() == _it.inf() && it.sup() == _it.sup()))
	{
		std::cerr.precision(20);
		std::cerr << "n" << it.inf() << ", " << it.sup() << std::endl;
		std::cerr << "_n " << _it.inf() << ", " << _it.sup() << std::endl;
		return false;
	}
	else
		return true;
}

TEST(LazyExactNumber, Operators)
{

	LazyNumber::Protector P;

	srand(0);
	for (size_t i = 0; i < 100000; i++)
	{
		int j = rand(), k = rand();
		if (j == 0 || k == 0)
			continue;

#if 1
		LazyNumber n1(j * M_PI);
		LazyNumber n2(k * M_PI);
		LazyNumber n3  = n1 + n2;
		LazyNumber n4  = n1 - n2;
		LazyNumber n5  = n1 * n2;
		LazyNumber n6  = n1 / n2;
		LazyNumber n7  = -(n4 + n5);
		LazyNumber n8  = -(n7 * n6);
		LazyNumber n9  = OMC::abs(n7 * n6);
		LazyNumber n10 = n9 / n5;
		LazyNumber n11 = n10 / n8;
		LazyNumber n12 = n10;
		if (n11 < n12)
			n12 = n11;

		n1.exact();
		n2.exact();
		n3.exact();
		n4.exact();
		n5.exact();
		n6.exact();
		n7.exact();
		n8.exact();
		n9.exact();
		n10.exact();
		n11.exact();
#endif

#if 1
		Lazy_exact_nt _n1(j * M_PI);
		Lazy_exact_nt _n2(k * M_PI);
		Lazy_exact_nt _n3  = _n1 + _n2;
		Lazy_exact_nt _n4  = _n1 - _n2;
		Lazy_exact_nt _n5  = _n1 * _n2;
		Lazy_exact_nt _n6  = _n1 / _n2;
		Lazy_exact_nt _n7  = -(_n4 + _n5);
		Lazy_exact_nt _n8  = -(_n7 * _n6);
		Lazy_exact_nt _n9  = CGAL::abs(_n7 * _n6);
		Lazy_exact_nt _n10 = _n9 / _n5;
		Lazy_exact_nt _n11 = _n10 / _n8;
		Lazy_exact_nt _n12 = _n10;
		if (_n11 < _n12)
			_n12 = _n11;

		_n1.exact();
		_n2.exact();
		_n3.exact();
		_n4.exact();
		_n5.exact();
		_n6.exact();
		_n7.exact();
		_n8.exact();
		_n9.exact();
		_n10.exact();
		_n11.exact();
#endif

#if 1
		ASSERT_TRUE(interval_number_equal(n1.approx(), _n1.approx()));
		ASSERT_TRUE(interval_number_equal(n2.approx(), _n2.approx()));
		ASSERT_TRUE(interval_number_equal(n3.approx(), _n3.approx()));
		ASSERT_TRUE(interval_number_equal(n4.approx(), _n4.approx()));
		ASSERT_TRUE(interval_number_equal(n5.approx(), _n5.approx()));
		ASSERT_TRUE(interval_number_equal(n6.approx(), _n6.approx()));
		ASSERT_TRUE(interval_number_equal(n7.approx(), _n7.approx()));
		ASSERT_TRUE(interval_number_equal(n8.approx(), _n8.approx()));
		ASSERT_TRUE(interval_number_equal(n9.approx(), _n9.approx()));
		ASSERT_TRUE(interval_number_equal(n10.approx(), _n10.approx()));
		ASSERT_TRUE(interval_number_equal(n11.approx(), _n11.approx()));
		ASSERT_TRUE(interval_number_equal(n12.approx(), _n12.approx()));

		ASSERT_EQ(n1.exact(), _n1.exact());
		ASSERT_EQ(n2.exact(), _n2.exact());
		ASSERT_EQ(n3.exact(), _n3.exact());
		ASSERT_EQ(n4.exact(), _n4.exact());
		ASSERT_EQ(n5.exact(), _n5.exact());
		ASSERT_EQ(n6.exact(), _n6.exact());
		ASSERT_EQ(n7.exact(), _n7.exact());
		ASSERT_EQ(n8.exact(), _n8.exact());
		ASSERT_EQ(n9.exact(), _n9.exact());
		ASSERT_EQ(n10.exact(), _n10.exact());
		ASSERT_EQ(n11.exact(), _n11.exact());
		ASSERT_EQ(n12.exact(), _n12.exact());
#endif
	}
}