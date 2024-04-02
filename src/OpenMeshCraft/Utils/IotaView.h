#pragma once

namespace OMC {
/**
 * @brief IotaView contains a sequential numbers from begin to (unbounded) end.
 * For example: from 1 to 10, it contains 1,2,3,4,5,6,7,8,9.
 * IotaView is a more simple class than std::ranges::iota_view in c++20.
 * In case when c++20 is unavailable, use IotaView to replace
 * std::ranges::iota_view.
 * @tparam T
 */
template <typename T>
class IotaView
{
	struct iterator
	{
	public:
		using iterator_category = std::random_access_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T *;
		using reference         = T &;

	public:
		// clang-format off
		iterator(T _num) : num(_num) { }
		T operator*() const { return num; }
		T operator*() { return num; }
		iterator &operator++() { ++num; return *this; }
		iterator &operator--() { --num; return *this; }
		iterator operator++(int) { iterator cpy(*this); num++; return cpy; }
		iterator operator--(int) { iterator cpy(*this); num--; return cpy; }
		bool operator==(const iterator &other) const { return num == other.num; }
		bool operator!=(const iterator &other) const { return num != other.num; }
		// clang-format on
	public:
		T num;
	};

	T begin_num;
	T end_num;

	IotaView()
	  : begin_num(0)
	  , end_num(0)
	{
	}

	IotaView(T _begin, T _end)
	  : begin_num(_begin)
	  , end_num(_end)
	{
	}

	iterator begin() { return iterator(begin_num); }
	iterator begin() const { return iterator(begin_num); }
	iterator cbegin() const { return iterator(begin_num); }
	iterator end() { return iterator(end_num); }
	iterator end() const { return iterator(end_num); }
	iterator cend() const { return iterator(end_num); }

	size_t size() const { return end_num - begin_num; }

	bool empty() const { return size() == 0; }
};
} // namespace OMC