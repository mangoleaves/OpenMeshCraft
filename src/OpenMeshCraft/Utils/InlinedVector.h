#pragma once

namespace OMC {
/// @brief It's similar to absl::inlined_vector, but more simple and more
/// friendly for debug.
/// @note Use this inlined vector to store simple types!!!
/// I don't know what will happen for complex types. :(
template <typename T, size_t N>
class InlinedVector
{
	// N is not zero and N is power of 2.
	static_assert(N != 0 && (N & (N - 1)) == 0);

	static_assert(std::is_trivially_copyable<T>::value);

public:
	using value_type = T;

public: /* Constructors *****************************************************/
	InlinedVector() noexcept
	  : cur_elem_size(0)
	{
	}

	InlinedVector(const InlinedVector &rhs) noexcept { operator=(rhs); }

	InlinedVector &operator=(const InlinedVector &rhs) noexcept
	{
		if (rhs.cur_elem_size <= N)
		{
			std::copy(rhs.inlined_data.begin(),
			          rhs.inlined_data.begin() + rhs.cur_elem_size,
			          inlined_data.begin());
			heap_data.clear();
		}
		else
		{
			inlined_data = rhs.inlined_data;
			heap_data    = rhs.heap_data;
		}
		cur_elem_size = rhs.cur_elem_size;
		return *this;
	}

	InlinedVector(InlinedVector &&rhs) noexcept { operator=(std::move(rhs)); }

	InlinedVector &operator=(InlinedVector &&rhs) noexcept
	{
		if (rhs.cur_elem_size <= N)
		{
			std::copy(rhs.inlined_data.begin(),
			          rhs.inlined_data.begin() + rhs.cur_elem_size,
			          inlined_data.begin());
			heap_data.clear();
		}
		else
		{
			inlined_data = rhs.inlined_data;
			heap_data    = std::move(rhs.heap_data);
		}
		cur_elem_size = rhs.cur_elem_size;
		return *this;
	}

public: /* Modifier **********************************************************/
	void push_back(const T &elem)
	{
		if (cur_elem_size < N)
			inlined_data[cur_elem_size] = elem;
		else
			heap_data.push_back(elem);
		cur_elem_size++;
	}
	void push_back(T &&elem)
	{
		if (cur_elem_size < N)
			inlined_data[cur_elem_size] = std::move(elem);
		else
			heap_data.push_back(elem);
		cur_elem_size++;
	}

	void pop_back()
	{
		if (cur_elem_size > N)
			heap_data.pop_back();
		cur_elem_size--;
	}

	void clear()
	{
		cur_elem_size = 0;
		heap_data.clear();
	}

	void resize(size_t s)
	{
		if (s > N) // need heap data
			heap_data.resize(s - N);
		else // s <= N, only inlined data is needed.
			heap_data.clear();
		cur_elem_size = s;
	}

	// void reserve(size_t s) {}

public: /* Access ************************************************************/
	const T &operator[](size_t idx) const
	{
		if (idx < N)
			return inlined_data[idx];
		else
			return heap_data[idx - N];
	}

	T &operator[](size_t idx)
	{
		if (idx < N)
			return inlined_data[idx];
		else
			return heap_data[idx - N];
	}

	T       &front() { return operator[](0); }
	const T &front() const { return operator[](0); }

	T       &back() { return operator[](cur_elem_size - 1); }
	const T &back() const { return operator[](cur_elem_size - 1); }

	size_t size() const { return (size_t)cur_elem_size; }

	bool empty() const { return cur_elem_size == 0; }

public: /* Iterator **********************************************************/
	      // These iterators are very simple and unsafe, use them carefully.
	      // clang-format off
	class iterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
	public:
		iterator(InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		iterator &operator++() { cur_idx++; return *this; }
		iterator &operator--() { cur_idx--; return *this; }
		iterator operator++(int) { iterator cpy(*this); cur_idx++; return cpy; }
		iterator operator--(int) { iterator cpy(*this); cur_idx--; return cpy; }
		iterator &operator+=(size_t step) { cur_idx += step; return *this; }
		iterator operator+(size_t step) const { iterator cpy(*this); cpy.cur_idx += step; return cpy; }
		iterator &operator-=(size_t step) { cur_idx -= step; return *this; }
		iterator operator-(size_t step) const { iterator cpy(*this); cpy.cur_idx -= step; return cpy; }
		size_t operator-(const iterator& other) const { return cur_idx - other.cur_idx; }
		T &operator*() const { return vec->operator[](cur_idx); }
		T &operator*() { return vec->operator[](cur_idx); }
		T *operator->() const { return &vec->operator[](cur_idx); }
		T *operator->() { return &vec->operator[](cur_idx); }
		bool operator<(const iterator &rhs) const { return cur_idx < rhs.cur_idx; }
		bool operator<=(const iterator &rhs) const { return cur_idx <= rhs.cur_idx; }
		bool operator>(const iterator &rhs) const { return cur_idx > rhs.cur_idx; }
		bool operator>=(const iterator &rhs) const { return cur_idx >= rhs.cur_idx; }
		bool operator==(const iterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const iterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		InlinedVector   *vec;
		difference_type  cur_idx;
	};

	iterator begin() {return iterator(this, 0);}
	iterator end() {return iterator(this, (typename iterator::difference_type)cur_elem_size);}

	class const_iterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = const T*;
		using reference         = const T&;
	public:
		const_iterator(const InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		const_iterator &operator++() { cur_idx++; return *this; }
		const_iterator &operator--() { cur_idx--; return *this; }
		const_iterator operator++(int) { const_iterator cpy(*this); cur_idx++; return cpy; }
		const_iterator operator--(int) { const_iterator cpy(*this); cur_idx--; return cpy; }
		const_iterator &operator+=(size_t step) { cur_idx += step; return *this; }
		const_iterator operator+(size_t step) { const_iterator cpy(*this); cpy.cur_idx += step; return cpy; }
		const_iterator &operator-=(size_t step) { cur_idx -= step; return *this; }
		const_iterator operator-(size_t step) { const_iterator cpy(*this); cpy.cur_idx -= step; return cpy; }
		size_t operator-(const const_iterator& other) const { return cur_idx - other.cur_idx; }
		const T &operator*() const { return vec->operator[](cur_idx); }
		const T &operator*() { return vec->operator[](cur_idx); }
		const T *operator->() const { return &vec->operator[](cur_idx); }
		const T *operator->() { return &vec->operator[](cur_idx); }
		bool operator<(const const_iterator &rhs) const { return cur_idx < rhs.cur_idx; }
		bool operator<=(const const_iterator &rhs) const { return cur_idx <= rhs.cur_idx; }
		bool operator>(const const_iterator &rhs) const { return cur_idx > rhs.cur_idx; }
		bool operator>=(const const_iterator &rhs) const { return cur_idx >= rhs.cur_idx; }
		bool operator==(const const_iterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const const_iterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		const InlinedVector   *vec;
		difference_type       cur_idx;
	};

	const_iterator cbegin() const {return const_iterator(this, 0);}
	const_iterator cend() const {return const_iterator(this, (typename const_iterator::difference_type)cur_elem_size);}
	const_iterator begin() const {return const_iterator(this, 0);}
	const_iterator end() const {return const_iterator(this, (typename const_iterator::difference_type)cur_elem_size);}

	class reverse_iterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
	public:
		reverse_iterator(InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		reverse_iterator &operator++() { cur_idx--; return *this; }
		reverse_iterator &operator--() { cur_idx++; return *this; }
		reverse_iterator operator++(int) { reverse_iterator cpy(*this); cur_idx--; return cpy; }
		reverse_iterator operator--(int) { reverse_iterator cpy(*this); cur_idx++; return cpy; }
		reverse_iterator &operator+=(size_t step) { cur_idx -= step; return *this; }
		reverse_iterator operator+(size_t step) const { reverse_iterator cpy(*this); cpy.cur_idx -= step; return cpy; }
		reverse_iterator &operator-=(size_t step) { cur_idx += step; return *this; }
		reverse_iterator operator-(size_t step) const { reverse_iterator cpy(*this); cpy.cur_idx += step; return cpy; }
		size_t operator-(const reverse_iterator& other) const { return -(cur_idx - other.cur_idx); }
		T &operator*() const { return vec->operator[](cur_idx); }
		T &operator*() { return vec->operator[](cur_idx); }
		T *operator->() const { return &vec->operator[](cur_idx); }
		T *operator->() { return &vec->operator[](cur_idx); }
		bool operator<(const reverse_iterator &rhs) const { return cur_idx > rhs.cur_idx; }
		bool operator<=(const reverse_iterator &rhs) const { return cur_idx >= rhs.cur_idx; }
		bool operator>(const reverse_iterator &rhs) const { return cur_idx < rhs.cur_idx; }
		bool operator>=(const reverse_iterator &rhs) const { return cur_idx <= rhs.cur_idx; }
		bool operator==(const reverse_iterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const reverse_iterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		InlinedVector   *vec;
		difference_type  cur_idx;
	};

	reverse_iterator rbegin() {return reverse_iterator(this, (typename reverse_iterator::difference_type)cur_elem_size - 1);}
	reverse_iterator rend() {return reverse_iterator(this, -1);}

	class const_reverse_iterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = const T*;
		using reference         = const T&;
	public:
		const_reverse_iterator(const InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		const_reverse_iterator &operator++() { cur_idx--; return *this; }
		const_reverse_iterator &operator--() { cur_idx++; return *this; }
		const_reverse_iterator operator++(int) { const_reverse_iterator cpy(*this); cur_idx--; return cpy; }
		const_reverse_iterator operator--(int) { const_reverse_iterator cpy(*this); cur_idx++; return cpy; }
		const_reverse_iterator &operator+=(size_t step) { cur_idx -= step; return *this; }
		const_reverse_iterator operator+(size_t step) const { const_reverse_iterator cpy(*this); cpy.cur_idx -= step; return cpy; }
		const_reverse_iterator &operator-=(size_t step) { cur_idx += step; return *this; }
		const_reverse_iterator operator-(size_t step) const { const_reverse_iterator cpy(*this); cpy.cur_idx += step; return cpy; }
		size_t operator-(const const_reverse_iterator& other) const { return -(cur_idx - other.cur_idx); }
		const T &operator*() const { return vec->operator[](cur_idx); }
		const T &operator*() { return vec->operator[](cur_idx); }
		const T *operator->() const { return &vec->operator[](cur_idx); }
		const T *operator->() { return &vec->operator[](cur_idx); }
		bool operator<(const const_reverse_iterator &rhs) const { return cur_idx > rhs.cur_idx; }
		bool operator<=(const const_reverse_iterator &rhs) const { return cur_idx >= rhs.cur_idx; }
		bool operator>(const const_reverse_iterator &rhs) const { return cur_idx < rhs.cur_idx; }
		bool operator>=(const const_reverse_iterator &rhs) const { return cur_idx <= rhs.cur_idx; }
		bool operator==(const const_reverse_iterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const const_reverse_iterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		const InlinedVector   *vec;
		difference_type  cur_idx;
	};

	const_reverse_iterator rbegin() const {return const_reverse_iterator(this, (typename reverse_iterator::difference_type)cur_elem_size - 1);}
	const_reverse_iterator rend() const {return const_reverse_iterator(this, -1);}
	const_reverse_iterator crbegin() const {return const_reverse_iterator(this, (typename const_reverse_iterator::difference_type)cur_elem_size - 1);}
	const_reverse_iterator crend() const {return const_reverse_iterator(this, -1);}
	      // clang-format on

public:
	std::array<T, N> inlined_data;
	std::vector<T>   heap_data;

	size_t cur_elem_size;
};
} // namespace OMC