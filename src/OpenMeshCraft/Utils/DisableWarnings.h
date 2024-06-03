#if defined(_MSC_VER)
	#pragma warning(push, 1)
	#pragma warning(disable : 4005 4702 4706 4717 4996)
#elif defined(__GNUC__)
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Warray-bounds"
	#pragma GCC diagnostic ignored "-Wstrict-aliasing"
	#pragma GCC diagnostic ignored "-Wignored-qualifiers"
	#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif