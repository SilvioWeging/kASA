# source: https://stackoverflow.com/a/29904501
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

project(kASA C CXX)

cmake_minimum_required(VERSION 2.8)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -march=native")

add_subdirectory(stxxl)
#add_subdirectory(zlib)
include_directories(${STXXL_INCLUDE_DIRS})
include_directories(zlib)

include(CheckCXXCompilerFlag)

# Thanks to 'Some programmer dude' on stackoverflow:

# Check for standard to use
check_cxx_compiler_flag(-std=c++17 HAVE_FLAG_STD_CXX17)
if(HAVE_FLAG_STD_CXX17)
    # Have -std=c++17, use it
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
else()
    check_cxx_compiler_flag(-std=c++11 HAVE_FLAG_STD_CXX11)
    if(HAVE_FLAG_STD_CXX11)
        # Have -std=c++11, use that
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    else()
        message(SEND_ERROR "No supported language standard detected. kASA needs at least C++11!")
    endif()
endif()



# on MacOSX with clang we need to use libc++ for C++11 headers
if (CMAKE_CXX_COMPILER MATCHES ".*clang[+][+]"
     OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   check_cxx_compiler_flag(-stdlib=libc++ CXX_HAS_STDLIB_LIBCXX)
   if(CXX_HAS_STDLIB_LIBCXX)
     set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
     #set(CMAKE_CXX_FLAGS "-D_GLIBCXX_USE_CXX11_ABI=0")
   else()
     message(SEND_ERROR "Compilation with clang requires libc++.")
   endif()
endif()


# if gcc or visual studio is used, openmp is not an issue
if(CMAKE_COMPILER_IS_GNUCXX OR MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
endif()

if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
    message("debug mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
endif()

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")


add_executable(kASA source/main.cpp zlib/gzstream.cpp)

target_link_libraries(kASA ${STXXL_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/zlib/zlibBuild/libz.a)
