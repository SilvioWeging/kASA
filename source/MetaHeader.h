/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2019 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#pragma once

#define kASA_VERSION 1.0

#include <iostream>
#include <cstdint>
#include <string>
#include <algorithm>
#include <tuple>
#include <bitset> // debug
#include <cmath>
#include <sstream>
#include <cstdio>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <limits>
#include <vector>
#include <iterator>
#include <numeric>
#include <list>
#include <limits.h>
#include <memory> //Smartpointer
#include <typeinfo>

//#include <emmintrin.h>
#include <omp.h>
#if _HAS_CXX17
	#if __has_include(<execution>)
		#include <execution>
	#endif
#endif

#include <stxxl/vector>
#include <stxxl/bits/stream/unique.h>
#include <stxxl/bits/containers/sorter.h>
#include <stxxl/algorithm>
//#include <stxxl/unordered_map>
//#include <stxxl/map>

using namespace std;

// Check windows
#if _WIN32 || _WIN64
#include <stxxl/bits/io/wincall_file.h>
	#if __has_include(<filesystem>)
		#include <filesystem>
	#else
		#include <experimental/filesystem>
	#endif
#include <zlib.h>
#include <gzstream.hpp>
#define _SILENCE_PARALLEL_ALGORITHMS_EXPERIMENTAL_WARNING

#pragma pack(push, 1)
struct packedPair { uint32_t first = 0; uint16_t second = 0; };
#pragma pack(pop)

#pragma pack(push, 1)
struct packedBigPair { 
	uint64_t first = 0; uint32_t second = 0; 
	packedBigPair& operator=(const tuple<uint64_t, uint32_t>& a) { 
		first = get<0>(a); 
		second = get<1>(a); 
		return *this; 
	} 
	
	bool operator==(const packedBigPair& a) const {
		return a.first == this->first && a.second == this->second; 
	} 
	
	bool operator<(const packedBigPair& b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}
};
#pragma pack(pop)

#pragma pack(push, 1)
struct packedBigPairTrie { 
	uint64_t second = 0; uint32_t first = 0; 

	packedBigPairTrie() { 
		second = 0; first = 0;
	} 
	
	packedBigPairTrie(const uint32_t& f, const uint64_t& s) : second(s), first(f) {}
};
#pragma pack(pop)

#if _WIN64
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif

// From zlib
#if defined(MSDOS) || defined(OS2) || defined(WIN32) || defined(__CYGWIN__)
#  include <fcntl.h>
#  include <io.h>
#  define SET_BINARY_MODE(file) setmode(fileno(file), O_BINARY)
#else
#  define SET_BINARY_MODE(file)
#endif

// Check GCC
#if __GNUC__
#include <stxxl/bits/io/syscall_file.h>
#include <dirent.h>
#include <gzstream.hpp>
struct __attribute__((packed)) packedPair { uint32_t first = 0; uint16_t second = 0; };

struct __attribute__((packed)) packedBigPair {
	uint64_t first = 0; uint32_t second = 0;
	packedBigPair& operator=(const tuple<uint64_t, uint32_t>& a) {
		first = get<0>(a);
		second = get<1>(a);
		return *this;
	}

	bool operator==(const packedBigPair& a) const {
		return a.first == this->first && a.second == this->second;
	}

	bool operator<(const packedBigPair& b) const {
		return (this->first < b.first || (!(b.first < this->first) && this->second < b.second));
	}
};

struct __attribute__((packed)) packedBigPairTrie {
	uint64_t second = 0; uint32_t first = 0;

	packedBigPairTrie() {
		second = 0;
		first = 0;
	}

	packedBigPairTrie(const uint32_t& f, const uint64_t& s) : second(s), first(f) {}
};

#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif


#if _WIN32 || _WIN64
typedef stxxl::wincall_file stxxlFile;
#define __PRETTY_FUNCTION__ __FUNCSIG__ 
#endif
#if __GNUC__
typedef stxxl::syscall_file stxxlFile;
#endif


typedef	stxxl::VECTOR_GENERATOR<packedBigPair, 4U, 4U, 2101248, stxxl::RC>::result contentVecType_32p; //2101248 is dividable by 4096(blocksize from stxxl) and 12(sizeof(packedBigPair)
typedef	stxxl::VECTOR_GENERATOR<packedBigPairTrie, 4U, 4U, 2101248, stxxl::RC>::result trieVector;
typedef stxxl::VECTOR_GENERATOR<packedPair, 4U, 4U, 2101248, stxxl::RC>::result index_t_p;

typedef	stxxl::VECTOR_GENERATOR<packedBigPair, 16U, 16U, 4096 * 12, stxxl::RC>::result contentVecType_32p_old;
//typedef	stxxl::VECTOR_GENERATOR<packedBigPairTrie, 16U, 16U, 4096 * 12, stxxl::RC>::result trieVector_old;
//typedef stxxl::VECTOR_GENERATOR<packedPair, 16U, 8U, 4096 * 6, stxxl::RC>::result index_t_p_old;
typedef stxxl::VECTOR_GENERATOR<packedPair, 4U, 4U, 2097156 * 4, stxxl::RC>::result index_t_p_old;
typedef stxxl::VECTOR_GENERATOR<packedBigPairTrie, 4U, 4U, 2097156, stxxl::RC>::result trieVector_old;

typedef uint32_t readIDType;