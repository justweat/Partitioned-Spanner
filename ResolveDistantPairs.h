//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H
#define PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H

#include "Utilities.h"
#include "SortPairs.h"
#include "ShortestPaths.h"

namespace spanners{

    vector<pair<size_t, size_t>> findCandidatePairs(const vector<Point> &points,
                                                    double t,
                                                    const vector<unique_ptr<Partition>>& leaves,
                                                    const Graph& leaderSpanner,
                                                    const set<pair<size_t, size_t>>& contiguousNeighborCheck){

        vector<pair<size_t, size_t>> validDistantPairs{};

        size_t n = leaves.size();

        if(n <= 1) return validDistantPairs;

        for(size_t i{}; i < n; ++i){
            for(size_t j{i + 1}; j < n; ++j){
                if(contiguousNeighborCheck.find(
                        make_pair(min(leaves[i]->partitionIndex, leaves[j]->partitionIndex),max(leaves[i]->partitionIndex, leaves[j]->partitionIndex)))
                   == contiguousNeighborCheck.end()){

                    vector<Point> uCorners{
                            Point{leaves[i]->bounding_box.xmin(), leaves[i]->bounding_box.ymin()},
                            Point{leaves[i]->bounding_box.xmin(), leaves[i]->bounding_box.ymax()},
                            Point{leaves[i]->bounding_box.xmax(), leaves[i]->bounding_box.ymin()},
                            Point{leaves[i]->bounding_box.xmax(), leaves[i]->bounding_box.ymax()}
                    };

                    vector<Point> vCorners{
                            Point{leaves[j]->bounding_box.xmin(), leaves[j]->bounding_box.ymin()},
                            Point{leaves[j]->bounding_box.xmin(), leaves[j]->bounding_box.ymax()},
                            Point{leaves[j]->bounding_box.xmax(), leaves[j]->bounding_box.ymin()},
                            Point{leaves[j]->bounding_box.xmax(), leaves[j]->bounding_box.ymax()}
                    };

                    bool validPair{};
                    for(const auto &u : uCorners){
                        for(const auto &v : vCorners){
                            if((getDistance(points[leaves[i]->leader], u) * t +

                               aStar(leaderSpanner.points, leaderSpanner.adjacencyMap, leaves[i]->leaderSpannerIndex, leaves[j]->leaderSpannerIndex) +

                               getDistance(points[leaves[j]->leader], v) * t) >

                               t * getDistance(u, v)){

                                validPair = true;
                                validDistantPairs.emplace_back(make_pair(leaves[i]->partitionIndex, leaves[j]->partitionIndex));
                                break;

                            }
                        }

                        if(validPair) break;

                    }
                }
            }
        }

        return validDistantPairs;

    }

    void mergeDistantPairs(const vector<Point>& points,
                           const unordered_map<size_t, vector<size_t>>& adjMap,
                           const vector<Point>& leaderPoints,
                           vector<unique_ptr<Partition>>& leaves,
                           const vector<pair<size_t, size_t>>& validDistantPairs,
                           const Graph& leaderSpanner,
                           double t,
                           size_t begin,
                           size_t end,
                           vector<Edge>& resolvedEdges,
                           mutex& resolvedEdgesLock){

        vector<Edge> threadEdges{};

        for(size_t i{begin}; i <= end; ++i){

            unique_ptr<Partition>& node_u = leaves[validDistantPairs[i].first];
            unique_ptr<Partition>& node_v = leaves[validDistantPairs[i].second];

            PointPairPQ pairs = sortPairs_DisjointCells(points, node_u->indices, node_v->indices);

            vector<tuple<size_t, size_t, number_t>> bridges;
            number_t initialBridge = aStar(points, adjMap, pairs.top().first.first, pairs.top().first.second);
            bridges.emplace_back(make_tuple(pairs.top().first.first, pairs.top().first.second, initialBridge));

            number_t leaderPathCost = aStar(leaderPoints, leaderSpanner.adjacencyMap, node_u->leaderSpannerIndex, node_v->leaderSpannerIndex);
            bridges.emplace_back(make_tuple(node_u->leader, node_v->leader, leaderPathCost));

            while(!pairs.empty()){

                auto [pts, dist] = pairs.top();
                size_t u{pts.first}, v{pts.second};
                pairs.pop();

                bool pathFound{};
                for(const auto &b : bridges){
                    auto [from, to, cost] = b;
                    if(
                            node_u->distances[node_u->localIndices.at(from)][node_u->localIndices.at(u)] +
                            cost +
                            node_v->distances[node_v->localIndices.at(to)][node_v->localIndices.at(v)] <
                            t * CGAL::sqrt(CGAL::squared_distance(points[u], points[v])))
                    {
                        pathFound = true;
                        break;
                    }
                }
                if(!pathFound){

                    bridges.emplace_back(make_tuple(u, v, CGAL::sqrt(CGAL::squared_distance(points[u], points[v]))));
                    threadEdges.emplace_back(Edge{u,v});

                }
            }
        }



        lock_guard<mutex> lk(resolvedEdgesLock);
        for(const auto& e : threadEdges){
            resolvedEdges.emplace_back(e);
        }

    }

    vector<Edge> resolveDistantPairs(const vector<Point> &points,
                                     const unordered_map<size_t, vector<size_t>>& adjMap,
                                     const vector<Point>& leaderPoints,
                                     vector<unique_ptr<Partition>>& leaves,
                                     double t,
                                     const Graph& leaderSpanner,
                                     const set<pair<size_t, size_t>>& contiguousNeighborCheck,
                                     unsigned int numOfThreads = 1){

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

        vector<pair<size_t, size_t>> validDistantPairs = findCandidatePairs(points,
                                                                            t,
                                                                            leaves,
                                                                            leaderSpanner,
                                                                            contiguousNeighborCheck);

        size_t n = validDistantPairs.size();

        size_t block = n / numOfThreads;
        if(block == 0){
            block = n;
        }

        vector<Edge> resolvedEdges{};
        mutex resolvedEdgesLock{};

        vector<thread> threads{};
        for(size_t i{}; i < n; i += block){
            threads.emplace_back(mergeDistantPairs,
                                 cref(points),
                                 cref(adjMap),
                                 cref(leaderPoints),
                                 ref(leaves),
                                 cref(validDistantPairs),
                                 cref(leaderSpanner),
                                 t,
                                 i,
                                 min(n - 1, i + block),
                                 ref(resolvedEdges),
                                 ref(resolvedEdgesLock));
        }

        for(auto& thread : threads){
            thread.join();
        }

        return resolvedEdges;

    }

}

#endif //PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H
