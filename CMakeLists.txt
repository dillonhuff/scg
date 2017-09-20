cmake_minimum_required(VERSION 2.6)

project(ecg)

SET( EXTRA_CXX_COMPILE_FLAGS "-std=c++11 -I./src -I./test -I/opt/local/include -O2")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_COMPILE_FLAGS}")

INCLUDE_DIRECTORIES(./src/)

SET(UTILS_CPPS ./src/utils/arena_allocator.cpp)

SET(UTILS_HEADERS ./src/utils/algorithm.h
		  ./src/utils/arena_allocator.h)

add_library(utils ${UTILS_CPPS} ${UTILS_HEADERS})

SET(GEOMETRY_CPPS ./src/geometry/point.cpp
		  ./src/geometry/surface.cpp
		  ./src/geometry/offset.cpp
		  ./src/geometry/polygon_3.cpp
		  ./src/geometry/extern_triangulate.cpp
	 ./src/geometry/arc.cpp
	 ./src/geometry/box.cpp
	 ./src/geometry/depth_field.cpp
	 #./src/geometry/extrusion.cpp
	 ./src/geometry/plane.cpp
	 ./src/geometry/line.cpp
	 ./src/geometry/matrix.cpp
	 ./src/geometry/polygon.cpp
	 ./src/geometry/polyline.cpp
	 ./src/geometry/spline_sampling.cpp
	 ./src/geometry/rotation.cpp
	 ./src/geometry/triangle.cpp
	 ./src/geometry/triangular_mesh.cpp
	 ./src/geometry/triangular_mesh_utils.cpp
	 ./src/geometry/trimesh.cpp
	 ./src/geometry/trimesh_types.cpp
	 ./src/system/parse_stl.cpp
	 ./src/geometry/homogeneous_transformation.cpp)

add_library(geometry ${GEOMETRY_CPPS})

target_link_libraries(geometry utils)

SET(GEOMETRY_TEST_FILES test/plane_tests.cpp
			test/rotation_tests.cpp
			test/matrix_tests.cpp
			test/surface_tests.cpp
			test/polygon_tests.cpp
			test/polyline_tests.cpp
			test/triangle_tests.cpp
			test/arc_tests.cpp
			test/ring_tests.cpp
			test/point_tests.cpp
			test/line_tests.cpp
			test/mesh_tests.cpp
			test/spline_tests.cpp)

add_executable(geometry-tests test/main_geometry.cpp ${GEOMETRY_TEST_FILES})
target_link_libraries(geometry-tests geometry utils)

