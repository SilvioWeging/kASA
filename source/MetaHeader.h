/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/
#pragma once

#define kASA_VERSION_MAJOR 1
#define kASA_VERSION_MINOR 4
#define kASA_VERSION_PATCH 4

#include <iostream>
#include <cstdint>
#include <string>
#include <algorithm>
#include <tuple>
#include <bitset> // debug
#include <cmath>
#include <sstream>
#include <cstdio>
#include <cstdlib>
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
#include <thread>
#include <mutex>
#include <queue>
#include <regex>

#include "utils/packedPairs.hpp"

//#include <immintrin.h> // AVX

//#include <omp.h>
#if _HAS_CXX17
	#if __has_include(<execution>)
		#include <execution>
	#endif
	//#if __has_include(<charconv>)
	//	#include <charconv>
	//#endif
#endif

#include <stxxl/vector>
#include <stxxl/bits/stream/unique.h>
#include <stxxl/bits/containers/sorter.h>
#include <stxxl/algorithm>
//#include <stxxl/unordered_map>
//#include <stxxl/map>

using namespace std;

// Check windows
#if (_WIN32 || _WIN64) && !__APPLE__
#include <stxxl/bits/io/wincall_file.h>
#include <Windows.h>
	#if __has_include(<filesystem>)
		#include <filesystem>
	#else
		#include <experimental/filesystem>
	#endif

	#if (__aarch64__ || _M_ARM)
		typename std::ifstream igzstream;
	#else
		#include <zlib.h>
		#include <gzstream.hpp>
	#endif
#define _SILENCE_PARALLEL_ALGORITHMS_EXPERIMENTAL_WARNING

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
#if __GNUC__ && !defined(__llvm__) && defined(_OPENMP)
#include <parallel/algorithm>
#endif

// Check GCC or Clang
#if (__GNUC__ || __clang__) && !_MSC_VER 
#include <sys/stat.h>
#include <unistd.h>
#include <stxxl/bits/io/syscall_file.h>
#include <dirent.h>
	#if (__aarch64__ || _M_ARM64)
		typename ifstream igzstream;
	#else
		#include <gzstream.hpp>
	#endif
#include <pthread.h>

#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32
#endif
#endif


#if (_WIN32 || _WIN64) && !__APPLE__
typedef stxxl::wincall_file stxxlFile;
#define __PRETTY_FUNCTION__ __FUNCSIG__ 
#endif
#if (__GNUC__ || defined(__llvm__)) && !_MSC_VER
typedef stxxl::syscall_file stxxlFile;
#endif

static bool _bShowLineDebugMode = false;
#define debugBarrier if (_bShowLineDebugMode ) { cerr << "File: " << __FILE__ << " Line: " << __LINE__ << endl; }


typedef	stxxl::VECTOR_GENERATOR<packedBigPair, 4U, 4U, 2101248, stxxl::RC>::result contentVecType_32p; //2101248 is dividable by 4096(blocksize from stxxl) and 12( sizeof(packedBigPair) )
typedef	stxxl::VECTOR_GENERATOR<packedLargePair, 4U, 4U, 2048000, stxxl::RC>::result contentVecType_128; // 2048000 = 4096*20*25
typedef	stxxl::VECTOR_GENERATOR<packedBigPairTrie, 4U, 4U, 2101248, stxxl::RC>::result trieVector;
typedef stxxl::VECTOR_GENERATOR<packedPair, 4U, 4U, 2101248, stxxl::RC>::result index_t_p;
typedef stxxl::VECTOR_GENERATOR<uint16_t, 4U, 4U, 2101248, stxxl::RC>::result taxaOnly;

typedef	stxxl::VECTOR_GENERATOR<packedBigPair, 16U, 16U, 4096 * 12, stxxl::RC>::result contentVecType_32p_old;
//typedef	stxxl::VECTOR_GENERATOR<packedBigPairTrie, 16U, 16U, 4096 * 12, stxxl::RC>::result trieVector_old;
//typedef stxxl::VECTOR_GENERATOR<packedPair, 16U, 8U, 4096 * 6, stxxl::RC>::result index_t_p_old;
typedef stxxl::VECTOR_GENERATOR<packedPair, 4U, 4U, 2097156 * 4, stxxl::RC>::result index_t_p_old;
typedef stxxl::VECTOR_GENERATOR<packedBigPairTrie, 4U, 4U, 2097156, stxxl::RC>::result trieVector_old;

typedef uint32_t readIDType;

constexpr auto HIGHESTPOSSIBLEK = 25;
constexpr uint64_t GIGABYTEASBYTES = 1024ull * 1024ull * 1024ull;

struct InputParameters {
	string cMode = "", sDBPathOut = "", sTempPath = "", sInput = "", contentFileIn = "", contentFile1 = "", contentFile2 = "", contentFileAfterUpdate = "", firstOldIndex = "", secondOldIndex = "", readToTaxaFile = "", tableFile = "", indexFile = "", delnodesFile = "", codonTable = "", sTaxonomyPath = "", sAccToTaxFiles = "", sTaxLevel = "", sStxxlMode = "", sCodonID = "1", sPairedEnd1 = "", sPairedEnd2 = "";
	bool bSpaced = false, bVerbose = false, bTranslated = false, bRAM = false, bUnique = false, bUnfunny = false, bSixFrames = false, bThreeFrames = false, bTaxIdsAsStrings = false, bCustomMemorySet = false, bVisualize = false, bHighKSetByUser = false, bOnlyOneFrame = false, bContinue = false, bCoverage = false;
	int32_t iNumOfThreads = 1, iHighestK = 12, iHigherK = 12, iLowerK = 7, iNumOfCall = 0, iNumOfBeasts = 3, iNumOfMask = 0;
	int64_t iMemorySizeAvail = 0;
	float fPercentageOfThrowAway = 0.f, threshold = 0.f;
	uint8_t iTrieDepth = 6, outputFormat = 1, shrinkStrategy = 1;//, iPrefixCheckMode = 0;
} GlobalInputParameters;