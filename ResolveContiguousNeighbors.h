//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H
#define PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H

#include "Utilities.h"
#include "QuadtreeNeighborFinder.h"
#include "SortPairs.h"

#include <thread>

namespace spanners{

    struct ContiguousNeighbors{
        vector<pair<size_t, size_t>> northNeighbors, eastNeighbors, neNeighbors, nwNeighbors;
        set<pair<size_t, size_t>> contiguousNeighborCheck{};
    };

    ContiguousNeighbors findContiguousNeighbors(const vector<QT_Node>& qtLeaves,
                                                const map<QT_Node, size_t>& leafIdentifier){

        ContiguousNeighbors contiguousNeighbors{};

        vector<Direction> contiguousMergeDirections{Direction::North, Direction::East, Direction::NorthEast, Direction::NorthWest};

        for(const auto& leaf : qtLeaves){
            for(auto direction : contiguousMergeDirections){
                QT_Node adjacent{};
                if(direction == Direction::North){

                    adjacent = leaf.adjacent_node(11);

                }else if(direction == Direction::East){

                    adjacent = leaf.adjacent_node(01);

                }else if(direction == Direction::NorthEast){

                    adjacent = leaf.adjacent_node(01);
                    if(!adjacent.is_null()){
                        adjacent = adjacent.adjacent_node(11);
                    }

                }else if(direction == Direction::NorthWest){

                    adjacent = leaf.adjacent_node(00);
                    if(!adjacent.is_null()){
                        adjacent = adjacent.adjacent_node(11);
                    }

                }

                if(!adjacent.is_null()){

                    vector<QT_Node> neighbors;
                    if(direction == Direction::North){

                        findLeafNeighbors_North(adjacent, neighbors);

                    }else if(direction == Direction::East){

                        findLeafNeighbors_East(adjacent, neighbors);

                    }else if(direction == Direction::NorthEast){

                        findLeafNeighbors_NE(adjacent, neighbors);

                    }else if(direction == Direction::NorthWest){

                        findLeafNeighbors_NW(adjacent, neighbors);

                    }

                    for(const auto& n : neighbors){

                        if(!n.is_null() && !n.empty()){

                            auto n1 = leafIdentifier.at(leaf);
                            auto n2 = leafIdentifier.at(n);

                            pair<size_t, size_t> neighborPair{min(n1, n2), max(n1, n2)};

                            contiguousNeighbors.contiguousNeighborCheck.insert(neighborPair);

                            if(direction == Direction::North){

                                contiguousNeighbors.northNeighbors.emplace_back(pair<size_t, size_t>{neighborPair});

                            }else if(direction == Direction::East){

                                contiguousNeighbors.eastNeighbors.emplace_back(pair<size_t, size_t>{neighborPair});

                            }else if(direction == Direction::NorthEast){

                                contiguousNeighbors.neNeighbors.emplace_back(pair<size_t, size_t>{neighborPair});

                            }else{

                                contiguousNeighbors.nwNeighbors.emplace_back(pair<size_t, size_t>{neighborPair});

                            }
                        }
                    }
                }
            }
        }

        return contiguousNeighbors;

    }

    void mergeContiguousNeighbors(const vector<Point>& points,
                                  vector<unique_ptr<Partition>>& leaves,
                                  const vector<pair<size_t, size_t>>& neighbors,
                                  number_t t,
                                  size_t begin,
                                  size_t end,
                                  Direction direction){

        for(size_t i{begin}; i < end; ++i){

            unique_ptr<Partition>& alphaInfo = leaves[neighbors[i].first];
            unique_ptr<Partition>& betaInfo = leaves[neighbors[i].second];

            vector<size_t>& alphaN = alphaInfo->indices;
            vector<size_t>& betaN = betaInfo->indices;

            PointPairPQ pairs = sortPairs_DisjointCells(points, alphaN, betaN);

            vector<pair<size_t, size_t>> bridges{};

            size_t counter{};

            if(alphaInfo->indices.size() == 1){
                alphaInfo->distances[0][0] = 0;
            }

            if(betaInfo->indices.size() == 1){
                betaInfo->distances[0][0] = 0;
            }

            while(!pairs.empty()){

                auto [pts, dist] = pairs.top();
                size_t u = pts.first, v = pts.second;
                pairs.pop();

                bool found_path{};
                for(const auto& b : bridges){

                    if(
                            alphaInfo->distances[alphaInfo->localIndices.at(u)][alphaInfo->localIndices.at(b.first)] +

                            getDistance(points[b.first], points[b.second]) +

                            betaInfo->distances[betaInfo->localIndices.at(v)][betaInfo->localIndices.at(b.second)] <

                            t * dist)
                    {

                        found_path = true;
                        break;

                    }
                }

                if(!found_path){
                    ++counter;
                    alphaInfo->contiguousResolutionEdges.emplace_back(Edge{u,v});
                    bridges.emplace_back(Edge{u, v});
                }
            }

//            cout << alphaInfo->indices.size() << "->" << betaInfo->indices.size() << " | " << counter << " | ";
//            if(direction == Direction::North){
//                cout << "North" << endl;
//            }else if(direction == Direction::East){
//                cout << "East" << endl;
//            }else if(direction == Direction::NorthEast){
//                cout << "NorthEast" << endl;
//            }else{
//                cout << "NorthWest" << endl;
//            }

        }
    }


    void resolveContiguousNeighbors(const vector<Point> &points,
                                    double t,
                                    vector<unique_ptr<Partition>>& leaves,
                                    const ContiguousNeighbors& neighbors,
                                    unsigned int numOfThreads = 1){

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

        if(!neighbors.northNeighbors.empty()){

            size_t n = neighbors.northNeighbors.size();
            size_t block = n / numOfThreads;
            if(block == 0){
                block = 1;
            }

            vector<thread> threads{};
            for(size_t i{}; i <= n; i += block){
                threads.emplace_back(mergeContiguousNeighbors, cref(points), ref(leaves), cref(neighbors.northNeighbors), t, i, min(n, i + block), Direction::North);
            }

            for(auto& thread : threads){
                thread.join();
            }

        }

        if(!neighbors.eastNeighbors.empty()){

            size_t n = neighbors.eastNeighbors.size();
            size_t block = ceil(n / numOfThreads);
            if(block == 0){
                block = n;
            }

            vector<thread> threads{};
            for(size_t i{}; i < n; i += block){
                threads.emplace_back(mergeContiguousNeighbors, cref(points), ref(leaves), cref(neighbors.eastNeighbors), t, i, min(n, i + block), Direction::East);
            }

            for(auto& thread : threads){
                thread.join();
            }

        }

        if(!neighbors.neNeighbors.empty()){

            size_t n = neighbors.neNeighbors.size();
            size_t block = ceil(n / numOfThreads);
            if(block == 0){
                block = n;
            }

            vector<thread> threads{};
            for(size_t i{}; i < n; i += block){
                threads.emplace_back(mergeContiguousNeighbors, cref(points), ref(leaves), cref(neighbors.neNeighbors), t, i, min(n, i + block), Direction::NorthEast);
            }

            for(auto& thread : threads){
                thread.join();
            }

        }

        if(!neighbors.nwNeighbors.empty()){

            size_t n = neighbors.nwNeighbors.size();
            size_t block = ceil(n / numOfThreads);
            if(block == 0){
                block = n;
            }

            vector<thread> threads{};
            for(size_t i{}; i < n; i += block){
                threads.emplace_back(mergeContiguousNeighbors, cref(points), ref(leaves), cref(neighbors.nwNeighbors), t, i, min(n, i + block), Direction::NorthWest);
            }

            for(auto& thread : threads){
                thread.join();
            }

        }

    }

}

#endif //PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H
