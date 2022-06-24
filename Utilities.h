#ifndef UNF_SPANNERS_UTILITIES_H
#define UNF_SPANNERS_UTILITIES_H

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>
#include <boost/functional/hash.hpp> // hashing pairs
#include <boost/heap/fibonacci_heap.hpp> // ordering

#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Vector_2.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/Orthtree/Node.h>
#include <CGAL/Quadtree.h>
#include <CGAL/Octree.h>


namespace spanners {

    using namespace std;

    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef K::Point_2 Point;
    typedef pair<size_t, size_t> Edge;
    typedef K::FT number_t;

    typedef CGAL::Quadtree<K, std::vector<Point>> Quadtree;
    typedef Quadtree::Node QT_Node;

    typedef priority_queue<pair<Edge,number_t>, vector<pair<Edge,number_t>>, function<bool(pair<Edge,number_t>&, pair<Edge,number_t>&)>> PointPairPQ;

    const number_t PI = 3.14159265358979323846264338327950288419716939937510;

    enum class PartitionedSpannerLevel {PartitionedSpanner, LeaderSpanner, AuxiliarySpanner};

    number_t getDistance(Point u, Point v){
        return CGAL::sqrt(CGAL::squared_distance(u, v));
    }

    struct Graph{
        vector<Point> points{};
        vector<Edge> edges{};
        unordered_map<size_t, vector<size_t>> adjacencyMap{};
    };

    struct Partition{
        Partition(const vector<size_t> &indices,
                  size_t index,
                  const vector<Point> &points){

            this->indices = indices;
            this->partitionIndex = index;
            this->points = points;

            for(const auto& i : indices){
                this->adjMap.insert(make_pair(i, vector<size_t>{}));
            }

        }
        Partition() = delete;
        size_t partitionIndex{};
        vector<size_t> indices{};
        vector<Point> points{};
        vector<Edge> greedyEdges{}, contiguousResolutionEdges{};
        unordered_map<size_t, vector<size_t>> adjMap;
        size_t leader{};
        size_t leaderSpannerIndex{};
        CGAL::Bbox_2 bounding_box{};
        vector<vector<number_t>> distances{};
        unordered_map<size_t, size_t> localIndices{};
    };

    struct CellInfo{
        vector<vector<number_t>> partitionDistanceMatrix{};
        unordered_map<size_t, size_t> partitionLocalIndices{};
        vector<Edge> partitionEdges{};
    };

    enum class Direction {North, East, South, West, NorthEast, NorthWest, SouthEast, SouthWest};



}

#endif // UNF_SPANNERS_UTILITIES_H
