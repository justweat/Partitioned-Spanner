//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_FINDCANDIDATEDISTANTPAIRS_H
#define PARTITIONEDSPANNER_FINDCANDIDATEDISTANTPAIRS_H

#include "Utilities.h"
#include "ShortestPaths.h"

namespace spanners{

    vector<pair<Aux_QT_Node*, Aux_QT_Node*>> findCandidatePairs(const vector<Point> &points,
                                                                double t,
                                                                const vector<Aux_QT_Node*>& leaves,
                                                                const Graph& leaderSpanner,
                                                                const set<pair<size_t, size_t>>& contig_neighbor_check){

        vector<pair<Aux_QT_Node*, Aux_QT_Node*>> validDistantPairs{};

        size_t n = leaves.size();

        if(n <= 1) return validDistantPairs;

        for(size_t i{}; i < n; ++i){
            for(size_t j{i + 1}; j < n; ++j){
                if(contig_neighbor_check.find(make_pair(
                        min(leaves[i]->index, leaves[j]->index),
                        max(leaves[i]->index, leaves[j]->index))) != contig_neighbor_check.end()) {
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

                    for(const auto &u : uCorners){
                        bool validPair{};
                        for(const auto &v : vCorners){
                            if(getDistance(points[leaves[i]->leader], u) +

                               aStar(points, leaderSpanner.adjacencyMap, leaves[i]->leaderActual, leaves[j]->leaderActual) +

                               getDistance(points[leaves[j]->leader], v) >

                               t * getDistance(u, v)){

                                validPair = true;
                                validDistantPairs.emplace_back(make_pair(leaves[i], leaves[j]));
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
}

#endif //PARTITIONEDSPANNER_FINDCANDIDATEDISTANTPAIRS_H
