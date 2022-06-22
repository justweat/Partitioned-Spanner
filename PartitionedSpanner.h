//
// Created by justin on 5/16/22.
//

#include "QuadtreeNeighborFinder.h"
#include "ResolveContiguousNeighbors.h"
#include "ResolveDistantPairs.h"
#include "InitializePartitions.h"

#ifndef SP_MER_17_PARTITIONEDSPANNER_H
#define SP_MER_17_PARTITIONEDSPANNER_H

namespace spanners{

    Graph partitionedSpanner(vector<Point> &points,
                             size_t cellSize,
                             double t,
                             bool leaderSpannerConstructor = false,
                             unsigned int numOfThreads = 1){

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

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

        vector<unique_ptr<Partition>> leaves{};
        unordered_map<size_t, Partition&> leafInfo{};
        map<QT_Node, size_t> leafIdentifier{};
        vector<QT_Node> qtLeaves{};

        size_t partitionIndex{};
        for(Quadtree::Node i : qt.traverse<CGAL::Orthtrees::Leaves_traversal>()){
            if(i.empty()) continue;
            vector<size_t> indices;
            vector<Point> pts;
            for(const auto &p : i){
                indices.emplace_back(PointToIndex.at(p));
                pts.emplace_back(p);
            }
            Partition auxNode{indices, partitionIndex, i, pts};
            leaves.push_back(make_unique<Partition>(auxNode));
            leafInfo.insert(pair<size_t, Partition&>{partitionIndex, auxNode});
            leafIdentifier.insert(pair<QT_Node, size_t>{i, partitionIndex});
            qtLeaves.emplace_back(i);
            ++partitionIndex;
        }

        thread iniPar{initializePartitions, ref(leaves), points, t, numOfThreads};
        ContiguousNeighbors contiguousNeighbors = findContiguousNeighbors(qtLeaves, leafIdentifier);
        iniPar.join();

        for(const auto& leaf : leaves){
            for(const auto& i : leaf->distances){
                for(const auto& j : i){
                    if(j == CGAL_IA_MAX_DOUBLE){
                        cout << "bad" << endl;
                    }
                }
            }
        }

        thread resConNei{resolveContiguousNeighbors, cref(points), t, ref(leaves), cref(contiguousNeighbors), numOfThreads};

        vector<Point> leaderPoints{};
        for(const auto& leaf : leaves){
            leaderPoints.push_back(IndexToPoint.at(leaf->leader));
        }

        Graph leaderSpanner = FG_GreedySpanner(leaderPoints, t);
//        if(!leaderSpannerConstructor){
//            leaderSpanner = partitionedSpanner(leaderPoints, cellSize, t, true, numOfThreads);
//        }

        unordered_map<size_t, vector<size_t>> adjMap;
        for(size_t i{}; i < n; ++i){
            adjMap.insert(make_pair(i, vector<size_t>{}));
        }

        resConNei.join();



        for(const auto& leaf : leaves){

            for(const auto& e : leaf->greedyEdges){
                adjMap.at(e.first).emplace_back(e.second);
                adjMap.at(e.second).emplace_back(e.first);
            }

            for(const auto& e : leaf->contiguousResolutionEdges){
                adjMap.at(e.first).emplace_back(e.second);
                adjMap.at(e.second).emplace_back(e.first);
            }

        }

        if(!leaderSpannerConstructor){
            for(const auto& e : leaderSpanner.edges){
                adjMap.at(e.first).emplace_back(e.second);
                adjMap.at(e.second).emplace_back(e.first);
            }
        }

        vector<Edge> distantEdges = resolveDistantPairs(points,
                                                        adjMap,
                                                        leaderPoints,
                                                        leaves,
                                                        t,
                                                        leaderSpanner,
                                                        contiguousNeighbors.contiguousNeighborCheck);

        vector<Edge> totalEdges{};

        size_t sumGre{}, sumRes{};

        for(const auto& leaf : leaves){

            for(const auto& e : leaf->greedyEdges){
                ++sumGre;
                totalEdges.emplace_back(e);
            }

            for(const auto& e : leaf->contiguousResolutionEdges){
                ++sumRes;
                totalEdges.emplace_back(e);
            }

            for(const auto& e : distantEdges){
                totalEdges.emplace_back(e);
            }

        }

        if(!leaderSpannerConstructor){
            for(const auto& e : leaderSpanner.edges){
                totalEdges.emplace_back(e);
            }
        }

        cout << "Greedy Edge Count: " << sumGre << endl;
        cout << "Contiguous Resolution Edge Count: " << sumRes << endl;
        cout << "Distant Resolution Edge Count: " << distantEdges.size() << endl;
        cout << "Leader Spanner Edge Count: " << leaderSpanner.edges.size() << endl;

        return Graph{points, totalEdges, adjMap};

    }

}

#endif //SP_MER_17_PARTITIONEDSPANNER_H
