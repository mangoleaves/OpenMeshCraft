#pragma once

namespace OMC {
/// @brief It's similar to absl::inlined_vector, but more simple.
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
	{ /* do nothing */
		OMC_THROW_NOT_IMPLEMENTED();
	}

	void reserve(size_t s)
	{ /* do nothing */
		OMC_THROW_NOT_IMPLEMENTED();
	}

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

	T       &back() { return operator[](cur_elem_size - 1); }
	const T &back() const { return operator[](cur_elem_size - 1); }

	size_t size() const { return (size_t)cur_elem_size; }

	bool empty() const { return cur_elem_size == 0; }

public: /* Iterator **********************************************************/
	      // These iterators are very simple and unsafe, use them carefully.
	      // clang-format off
	class NoSkipIterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = T*;
		using reference         = T&;
	public:
		NoSkipIterator(InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		NoSkipIterator &operator++() { cur_idx++; return *this; }
		NoSkipIterator &operator--() { cur_idx--; return *this; }
		NoSkipIterator operator++(int) { NoSkipIterator cpy(*this); cur_idx++; return cpy; }
		NoSkipIterator operator--(int) { NoSkipIterator cpy(*this); cur_idx--; return cpy; }
		T &operator*() const { return vec->operator[](cur_idx); }
		T &operator*() { return vec->operator[](cur_idx); }
		T *operator->() const { return &vec->operator[](cur_idx); }
		T *operator->() { return &vec->operator[](cur_idx); }
		bool operator==(const NoSkipIterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const NoSkipIterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		InlinedVector   *vec;
		difference_type  cur_idx;
	};
	using iterator = NoSkipIterator;

	NoSkipIterator begin() {return NoSkipIterator(this, 0);}
	NoSkipIterator end() {return NoSkipIterator(this, (typename NoSkipIterator::difference_type)cur_elem_size);}

	class NoSkipConstIterator
	{
	public:
		using iterator_category = std::bidirectional_iterator_tag;
		using difference_type   = std::ptrdiff_t;
		using value_type        = T;
		using pointer           = const T*;
		using reference         = const T&;
	public:
		NoSkipConstIterator(const InlinedVector *_vec, difference_type _cur_idx)
		  : vec(_vec) , cur_idx(_cur_idx) {}
		NoSkipConstIterator &operator++() { cur_idx++; return *this; }
		NoSkipConstIterator &operator--() { cur_idx--; return *this; }
		NoSkipConstIterator operator++(int) { NoSkipConstIterator cpy(*this); cur_idx++; return cpy; }
		NoSkipConstIterator operator--(int) { NoSkipConstIterator cpy(*this); cur_idx--; return cpy; }
		const T &operator*() const { return vec->operator[](cur_idx); }
		const T &operator*() { return vec->operator[](cur_idx); }
		const T *operator->() const { return &vec->operator[](cur_idx); }
		const T *operator->() { return &vec->operator[](cur_idx); }
		bool operator==(const NoSkipConstIterator &rhs) const { return cur_idx == rhs.cur_idx; }
		bool operator!=(const NoSkipConstIterator &rhs) const { return cur_idx != rhs.cur_idx; }
	public:
		const InlinedVector   *vec;
		difference_type       cur_idx;
	};
	using const_iterator = NoSkipConstIterator;

	NoSkipConstIterator cbegin() const {return NoSkipConstIterator(this, 0);}
	NoSkipConstIterator cend() const {return NoSkipConstIterator(this, (typename NoSkipConstIterator::difference_type)cur_elem_size);}
	NoSkipConstIterator begin() const {return NoSkipConstIterator(this, 0);}
	NoSkipConstIterator end() const {return NoSkipConstIterator(this, (typename NoSkipConstIterator::difference_type)cur_elem_size);}
	      // clang-format on

public:
	std::array<T, N> inlined_data;
	std::vector<T>   heap_data;

	size_t cur_elem_size;
};
} // namespace OMC