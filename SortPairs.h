//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_SORTPAIRS_H
#define SP_MER_17_SORTPAIRS_H

#include "Utilities.h"

namespace spanners{

    function<bool(const pair<Edge, number_t>&, const pair<Edge, number_t>&)> PointPairPQ_comparator =
            [](const pair<Edge, number_t> &o1, const pair<Edge, number_t> &o2)->bool{
                return o2.second < o1.second;
            };

    /*
     * Sorts all n choose 2 points into non-decreasing order
     */
    PointPairPQ sortPairs(const vector<Point> &points, const vector<size_t> &indices){

        vector<pair<Edge, number_t>> pairs;
        size_t n = indices.size();

        for(size_t i{}; i < n; ++i){
            for(size_t j {i + 1}; j < n; ++j){
                pairs.emplace_back(make_pair(Edge{indices[i], indices[j]}, CGAL::sqrt(CGAL::squared_distance(points[indices[i]], points[indices[j]]))));
            }
        }
        return PointPairPQ(PointPairPQ_comparator, pairs);
    }

    /*
     * Sorts all points from node1 to node2 in non-decreasing order
     */
    PointPairPQ sortPairs_DisjointCells(const vector<Point> &points, const vector<size_t> &node1, const vector<size_t> &node2){

        vector<pair<Edge, number_t>> pairs;
        size_t n = node1.size();

        for(size_t i{}; i < n; ++i){
            for(const auto &j : node2){
                pairs.emplace_back(make_pair(Edge{node1[i], j}, CGAL::sqrt(CGAL::squared_distance(points[node1[i]], points[j]))));
            }
        }
        return PointPairPQ(PointPairPQ_comparator, pairs);
    }

}

#endif //SP_MER_17_SORTPAIRS_H
