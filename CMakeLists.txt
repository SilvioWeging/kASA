project(kASA)


cmake_minimum_required(VERSION 2.8)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(SEND_ERROR "In-source builds are not allowed.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

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
    check_cxx_compiler_flag(-std=c++14 HAVE_FLAG_STD_CXX14)
    if(HAVE_FLAG_STD_CXX14)
        # Have -std=c++14, use it
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    else()
        message(SEND_ERROR "No supported language standard detected. kASA needs at least C++14!")
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")


add_executable(kASA source/main.cpp zlib/gzstream.cpp)

if(APPLE)
	target_link_libraries(kASA ${STXXL_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/zlib/zlibBuild/libz.1.2.11.dylib)
else()
	target_link_libraries(kASA ${STXXL_LIBRARIES} ${CMAKE_CURRENT_SOURCE_DIR}/zlib/zlibBuild/libz.so.1.2.11)
endif(APPLE)