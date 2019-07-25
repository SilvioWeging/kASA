// -*- mode: c++ -*-
/***************************************************************************
 *  include/stxxl/bits/config.h.in
 *
 *  Template file processed by cmake to set all define switches for this build
 *  according to the cmake build options.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
 *  Modified      2019 Silvio Weging <silvio.weging@gmail.com>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONFIG_HEADER
#define STXXL_CONFIG_HEADER

// the STXXL library version variables
#define STXXL_VERSION_MAJOR 1
#define STXXL_VERSION_MINOR 4
#define STXXL_VERSION_PATCH 1
#define STXXL_VERSION_STRING "1.4.1"
#define STXXL_VERSION_PHASE "prerelease/"

// if this is a git repository, add the refspec and commit sha
/* #undef STXXL_VERSION_GIT_REFSPEC */
/* #undef STXXL_VERSION_GIT_SHA1 */

/* #undef STXXL_DIRECT_IO_OFF */
// default: 0/1 (platform dependent)
// cmake:   detection of platform and flag
// effect:  disables use of O_DIRECT flag on unsupported platforms

/* #undef STXXL_HAVE_MMAP_FILE */
// default: 0/1 (platform dependent)
// used in: io/mmap_file.h/cpp
// effect:  enables/disables memory mapped file implementation

/* #undef STXXL_HAVE_LINUXAIO_FILE */
// default: 0/1 (platform dependent)
// used in: io/linuxaio_file.h/cpp
// effect:  enables/disables Linux AIO file implementation

/* #undef STXXL_POSIX_THREADS */
// default: off
// cmake:   detection of pthreads by cmake
// effect:  uses POSIX thread and mutexes on linux

#define STXXL_STD_THREADS 1
// default: off
// cmake:   detection of <thread> and <mutex>
// effect:  uses std thread and mutex on windows or (forced on linux)

#define STXXL_WINDOWS 1
// default: off
// cmake:   detection of ms windows platform (32- or 64-bit)
// effect:  enables windows-specific api calls (mingw or msvc)

#define STXXL_MSVC 1910
// default: off
// cmake:   detection of ms visual c++ via CMake (contains version number)
// effect:  enables msvc-specific headers and macros

#define STXXL_HAVE_CXX11_RANGE_FOR_LOOP 1
// default: off
// run-time: detection C++11 support for "for (auto i : obj) { }"
// effect:  enables some C++11 construct (currently only allowed in examples)

/* #undef STXXL_HAVE_SYNC_ADD_AND_FETCH */
// default: off
// cmake:   detection of __sync_add_and_fetch() intrinsic
// effect:  enables use of atomics in counting_ptr

/* #undef STXXL_PARALLEL_MODE_EXPLICIT */
// default: off
// cmake:   -DUSE_GNU_PARALLEL=ON
// effect:  explicitly enables use of __gnu_parallel algorithms

/* #undef STXXL_BOOST_CONFIG */
/* #undef STXXL_BOOST_FILESYSTEM */
/* #undef STXXL_BOOST_RANDOM */
/* #undef STXXL_BOOST_THREADS */
/* #undef STXXL_BOOST_TIMESTAMP */
// default: off
// cmake:   -DUSE_BOOST=ON
// effect:  enables use of boost libraries in different parts of STXXL.

#if STXXL_BOOST_CONFIG
  #include <boost/config.hpp>
#endif

#define STXXL_STD_RANDOM 1
// default: off
// cmake:   detection of <random>
// effect:  uses std random generator on windows or (forced on linux)

/* #undef STXXL_HAVE_MALLINFO_PROTO */
// default: off
// cmake:   detection of mallinfo() function in <malloc.h>
// effect:  used by stxxl_tool/mallinfo for malloc stats

/* #undef STXXL_HAVE_MLOCK_PROTO */
// default: off
// cmake:   detection of mlock() function in <sys/mman.h>
// effect:  used by stxxl_tool/mlock for locking physical pages

/* #undef STXXL_WITH_VALGRIND */
// default: off
// cmake:   option USE_VALGRIND=ON
// effect:  run all tests with valgrind and pre-initialize some memory buffers

#endif // !STXXL_CONFIG_HEADER
