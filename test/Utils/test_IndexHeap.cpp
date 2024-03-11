#ifndef OMC_ENABLE_EXPENSIVE_ASSERT
#define OMC_ENABLE_EXPENSIVE_ASSERT
#endif

#include "OpenMeshCraft/Utils/IndexHeap.h"

#include "test_utils.h"

using namespace OMC;

TEST(Index, Heap)
{
	class Test_IndexDenseHeap : public IndexDenseHeap<int>
	{
	public:
		void print()
		{
			for (index_t i = 0; i < m_pos_in_heap.size(); i++)
			{
				if (is_valid_idx(m_pos_in_heap[i]))
					std::cout << m_pos_in_heap[i] << " ";
				else
					std::cout << "X ";
			}
			std::cout << std::endl;
			for (index_t i = 1; i < m_heap.size(); i++)
			{
				std::cout << std::format("({}, {}) ", m_heap[i].first,
				                         m_heap[i].second);
			}
			std::cout << std::endl;
		}
	};
	Test_IndexDenseHeap heap;
	heap.push(1, 3);
	EXPECT_ANY_THROW(heap.push</*AllowUpdate*/ false>(1, 2));
	heap.push</*AllowUpdate*/ true>(1, 2);
	heap.push(2, 4);
	heap.push(4, 1);
	heap.push(3, 3);
	EXPECT_ANY_THROW(heap.push</*AllowUpdate*/ false>(3, 5));
	heap.push</*AllowUpdate*/ true>(3, 5);
	heap.push(6, 9);
	heap.push(8, 6);
	heap.push(10, 0);

	ASSERT_EQ(heap.exist(InvalidIndex), false);
	ASSERT_EQ(heap.exist(10), true);
	ASSERT_EQ(heap.exist(5), false);
	ASSERT_EQ(heap.exist(0), false);
	ASSERT_EQ(heap.exist(1), true);

	ASSERT_EQ(heap.exist(8), true);
	heap.remove(8);
	ASSERT_EQ(heap.exist(8), false);

	ASSERT_EQ(heap.top().first, 10);
	heap.remove(10);
	ASSERT_EQ(heap.exist(10), false);

	ASSERT_EQ(heap.top().first, 4);
	ASSERT_EQ(heap.top().second, 1);
	heap.pop();
	ASSERT_EQ(heap.top().first, 1);
	ASSERT_EQ(heap.top().second, 2);
	heap.pop();
	ASSERT_EQ(heap.top().first, 2);
	ASSERT_EQ(heap.top().second, 4);
	heap.pop();
	ASSERT_EQ(heap.top().first, 3);
	ASSERT_EQ(heap.top().second, 5);
	heap.pop();
	ASSERT_EQ(heap.top().first, 6);
	ASSERT_EQ(heap.top().second, 9);
	heap.pop();
	ASSERT_EQ(heap.heap_size(), 0);
}