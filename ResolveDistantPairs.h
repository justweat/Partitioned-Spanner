//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H
#define PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H

#include "Utilities.h"
#include "SortPairs.h"
#include "ShortestPaths.h"

namespace spanners{

    /*
     * Candidate distant pair defined as:
     * From Partition_i to Partition_j, if:
     * distance from centroid of Partition_i to any bounding box corner of Partition_i * t +
     * leader spanner cost from Partition_i to Partition_j +
     * distance from centroid of Partition_j to any bounding box corner of Partition_j * t >
     * d(bounding box corner from Partition_i, bounding box corner from Partition_j)
     * consider this a candidate pair to be resolved
     *
     * Params:
     * points: entire point set
     * t: spanner invariant
     * leaves: unique ptrs to each Partition
     * leaderSpanner: centroid spanner used to ensure connected graph
     * contiguousNeighborCheck: determines whether this pair has already been resolved contiguously
     *
     */
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

    /*
     * Merges distant pairs using the Bridging heurstic described in PartitionedSpanners.h
     *
     * Params:
     * points: entire point set
     * adjMap: adjacency list of entire graph
     * leaderPoints: centroids from partitions
     * leaves: unique pts to each partition
     * validDistantPairs: candidate pairs from findCandidatePairs()
     * leaderSpanner: centroid spanner
     * t: spanner invariant
     * begin/end: indices into validDistantPairs [begin, end]
     * resolvedEdges: all edges used to resolve distant pairs between all threads
     * resolvedEdgesLock: mutex for resolvedEdges
     *
     * Returns:
     * void
     */
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

    /*
     * Launch threads to concurrently resolve distant pairs
     *
     * Params:
     * points: entire point set
     * adjMap: adjacency list of entire graph
     * leaderPoints: centroids from partitions
     * leaves: unique pts to each partition
     * t: spanner invariant
     * leaderSpanner: centroid spanner
     * contiguousNeighborCheck:set of partition indices that were involved in contiguous neighbor resolution
     * numOfThreads: max number of threads for concurrent resolution
     *
     * Returns:
     * all edges needed to resolve distant pairs
     *
     */
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
        for(size_t i{}; i < n; i += block + 1){
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
