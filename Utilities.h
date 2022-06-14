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

    function<bool(const pair<Edge, number_t>&, const pair<Edge, number_t>&)> PointPairPQ_comparator =
            [](const pair<Edge, number_t> &placed, const pair<Edge, number_t> &searching)->bool{
                return searching.second < placed.second;
            };

    const number_t PI = 3.14159265358979323846264338327950288419716939937510;

    number_t getDistance(Point u, Point v){
        return CGAL::sqrt(CGAL::squared_distance(u, v));
    }

    struct Graph{
        vector<Point> points{};
        vector<Edge> edges{};
        unordered_map<size_t, vector<size_t>> adjacencyMap{};
    };

    struct Aux_QT_Node{
        Aux_QT_Node(const vector<size_t> &indices, size_t index,
                    CGAL::Quadtree<spanners::K,
                    std::vector<spanners::Point>>::Node node,
                    const vector<Point> &points){
            this->indices = indices;
            this->index = index;
            this->original_node = &node;
            this->points = points;
        }
        Aux_QT_Node() = default;
        size_t index{};
        vector<size_t> indices{};
        vector<Point> points{};
        vector<Edge> edges{};
        size_t leader{};
        size_t leaderActual{};
        CGAL::Bbox_2 bounding_box{};
        CGAL::Quadtree<spanners::K, std::vector<spanners::Point>>::Node* original_node{};
        vector<vector<number_t>> distances{};
        unordered_map<size_t, size_t> local_indices{};
    };

    struct CellInfo{
        vector<vector<number_t>> distances;
        unordered_map<size_t, size_t> indices;
        vector<Edge> edges;
    };

    enum class Direction {North, East, South, West, NorthEast, NorthWest, SouthEast, SouthWest};



}

#endif // UNF_SPANNERS_UTILITIES_H
