//
// Created by justin on 5/16/22.
//

#include "QuadtreeNeighborFinder.h"
#include "ResolveContiguousNeighbors.h"
#include "FindCentroids.h"
#include "FindCandidateDistantPairs.h"
#include "ResolveValidDistantPairs.h"
#include "GreedySpanners.h"

#ifndef SP_MER_17_PARTITIONEDSPANNER_H
#define SP_MER_17_PARTITIONEDSPANNER_H

namespace spanners{

    Graph partitionedSpanner(vector<Point> &points, size_t cellSize, double t, bool leader = false){

        //prevents pass by const ref
        //TODO: investigate whether altering CGAL implementation is feasible
        Quadtree qt(points);
        qt.refine(numeric_limits<size_t>::max(), cellSize);

        size_t n = points.size();

        unordered_map<Point, size_t> PointToIndex{};
        unordered_map<size_t, Point> IndexToPoint{};

        for(size_t i{}; i < n; ++i){
            PointToIndex.insert(pair<Point, size_t>{points[i], i});
            IndexToPoint.insert(pair<size_t, Point>{i, points[i]});
        }

        unordered_map<size_t, Aux_QT_Node> leaf_info{};
        map<QT_Node, size_t> leaf_identifier{};
        vector<QT_Node> qt_leaves{};

        size_t leaf_index{}, num_of_leaves{};
        for(Quadtree::Node i : qt.traverse<CGAL::Orthtrees::Leaves_traversal>()){
            if(i.empty()) continue;
            qt_leaves.emplace_back(i);
            vector<size_t> indices;
            vector<Point> pts;
            for(const auto &p : i){
                indices.emplace_back(PointToIndex.at(p));
                pts.emplace_back(p);
            }
            leaf_info.insert(pair<size_t, Aux_QT_Node>{leaf_index, Aux_QT_Node{indices, leaf_index, i, pts}});
            leaf_identifier.insert(pair<QT_Node, size_t>{i, leaf_index});
            ++leaf_index; ++num_of_leaves;
        }

        unordered_map<size_t, vector<size_t>> adjMap;
        for(size_t i{}; i < n; ++i){
            adjMap.insert(make_pair(i, vector<size_t>{}));
        }

        vector<Edge> totalEdges{};

        for(auto &i : leaf_info){
            CellInfo cur = FG_GreedySpannerPartitionedCell(points, i.second.indices, adjMap, t);
            i.second.edges = cur.edges;
            i.second.distances = cur.distances;
            i.second.local_indices = cur.indices;
            for(const auto &e : i.second.edges){
                totalEdges.emplace_back(e);
            }
        }

        set<pair<size_t, size_t>> contiguous_neighbor_check;
        vector<Edge> resolved_edges{};

        vector<Direction> directions{Direction::North, Direction::East, Direction::NorthEast, Direction::NorthWest};

        for(auto &node : qt_leaves){
            for(auto dir : directions){
                vector<Edge> temp = resolveContiguousNeighbors(points,
                                                               t,
                                                               node,
                                                               leaf_info,
                                                               leaf_identifier,
                                                               dir,
                                                               contiguous_neighbor_check);
                for(const auto &e : temp){
                    adjMap.at(e.first).emplace_back(e.second);
                    adjMap.at(e.second).emplace_back(e.first);
                    totalEdges.emplace_back(e);
                }
            }
        }

        for(auto &i : leaf_info){
            i.second.bounding_box = CGAL::bbox_2(i.second.points.begin(), i.second.points.end());
        }

        vector<Point> leaders = findCentroid(points, leaf_info, IndexToPoint);

        Graph leaderSpanner;
        if(!leader){
            leaderSpanner = partitionedSpanner(leaders, cellSize, t, true);
            for(const auto &e : leaderSpanner.edges){
                totalEdges.emplace_back(e);
            }
        }

        vector<Aux_QT_Node*> leaves;
        for(auto &i : leaf_info){
            leaves.emplace_back(&i.second);
        }

        vector<pair<Aux_QT_Node*, Aux_QT_Node*>> validDistantPairs = findCandidatePairs(points,
                                                                                        t,
                                                                                        leaves,
                                                                                        leaderSpanner,
                                                                                        contiguous_neighbor_check);

        for(const auto &i : validDistantPairs){
            vector<Edge> temp = resolveValidDistantPairs(points,
                                                         leaders,
                                                         i.first,
                                                         i.second,
                                                         t,
                                                         leaderSpanner);

            for(const auto &e : temp){
                adjMap.at(e.first).emplace_back(e.second);
                adjMap.at(e.second).emplace_back(e.first);
                totalEdges.emplace_back(e);
            }

        }

        return Graph{points, totalEdges, adjMap};

    }

}

#endif //SP_MER_17_PARTITIONEDSPANNER_H
