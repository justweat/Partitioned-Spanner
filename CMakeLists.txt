cmake_minimum_required(VERSION 3.21)
project(Fast_Sparse_Spanner)

set(CMAKE_CXX_STANDARD 17)

add_executable(Fast_Sparse_Spanner main.cpp PointGenerators.h GraphPrinter.h LatexPrinter.h Utilities.h GreedySpanners.h SortPairs.h ShortestPaths.h Testing.h PointGenOptions.h QuadtreeNeighborFinder.h PartitionedSpanner.h ResolveContiguousNeighbors.h FindCentroids.h FindCandidateDistantPairs.h ResolveValidDistantPairs.h)

set(CMAKE_CXX_FLAGS "-O0 -fopenmp -fexceptions -frounding-math -Wall -Wextra -Wno-unused")

target_link_libraries(Fast_Sparse_Spanner gmp)
target_link_libraries(Fast_Sparse_Spanner mpfr)
target_link_libraries(Fast_Sparse_Spanner boost_system)
target_link_libraries(Fast_Sparse_Spanner boost_thread)
target_link_libraries(Fast_Sparse_Spanner -I/usr/local/include/CGAL)