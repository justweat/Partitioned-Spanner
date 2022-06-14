//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H
#define PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H

#include "Utilities.h"
#include "SortPairs.h"
#include "ShortestPaths.h"

namespace spanners{

    vector<Edge> resolveValidDistantPairs(const vector<Point> &points,
                                          const vector<Point>& leaders,
                                          Aux_QT_Node* node_u,
                                          Aux_QT_Node* node_v,
                                          double t,
                                          const Graph& leaderSpanner){

        vector<Edge> resolvedEdges{};

        PointPairPQ pairs = sortPairs_DisjointCells(points, node_u->indices, node_v->indices);

        vector<tuple<size_t, size_t, number_t>> bridges;
        number_t leaderPathCost = aStar(leaders, leaderSpanner.adjacencyMap, node_u->leaderActual, node_v->leaderActual);
        bridges.emplace_back(make_tuple(node_u->leader, node_v->leader, leaderPathCost));

        while(!pairs.empty()){

            auto [pts, dist] = pairs.top();
            size_t u{pts.first}, v{pts.second};
            pairs.pop();

            bool pathFound{};
            for(const auto &b : bridges){
                auto [from, to, cost] = b;
                if(
                        node_u->distances[node_u->local_indices.at(from)][node_u->local_indices.at(u)] +
                        cost +
                        node_v->distances[node_v->local_indices.at(to)][node_v->local_indices.at(v)] <
                        t * CGAL::sqrt(CGAL::squared_distance(points[u], points[v])))
                {
                    pathFound = true;
                    break;
                }
            }
            if(!pathFound){

                bridges.emplace_back(make_tuple(u, v, CGAL::sqrt(CGAL::squared_distance(points[u], points[v]))));
                resolvedEdges.emplace_back(Edge{u,v});

            }
        }

        return resolvedEdges;

    }

}

#endif //PARTITIONEDSPANNER_RESOLVEVALIDDISTANTPAIRS_H
