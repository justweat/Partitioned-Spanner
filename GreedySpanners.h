//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_GREEDYSPANNERS_H
#define SP_MER_17_GREEDYSPANNERS_H

#include "Utilities.h"
#include "ShortestPaths.h"
#include "SortPairs.h"


namespace spanners{

    /*
     * Constructs a geometric spanner where all |V| choose 2 pairs (u,v)
     * within the graph contain the invariant that their
     * graphPath(u,v) <= t * d(u,v)
     *
     * Specialized for partitioned cells
     *
     * Params:
     * points: entire point set
     * indices: original index for each point within the point set
     * adjMap: adjacency list for this partition
     * t: spanner invariant
     *
     * Returns:
     * CellInfo consisting of distance matrix for this partition, indices used to index into distance matrix, and edges
     */
    CellInfo FG_GreedySpanner_PartitionedCell(const vector<Point> &points,
                                              const vector<size_t> &indices,
                                              unordered_map<size_t, vector<size_t>> &adjMap,
                                              double t){

        vector<Edge> edges{};
        size_t n = indices.size();

        unordered_map<size_t, size_t> local_indices{};
        for(size_t i{}; i < n; ++i){
            local_indices.insert(make_pair(indices[i], i));
        }

        PointPairPQ pq = sortPairs(points, indices);

        vector<vector<number_t>> distances(n, vector<number_t>(n, CGAL_IA_MAX_DOUBLE));

        for(size_t i{}; i < n; ++i){
            distances[i][i] = 0;
        }

        while(!pq.empty()){

            auto [pts, euc_u_to_v] = pq.top();
            auto u = pts.first, v = pts.second;
//            euc_u_to_v = CGAL::sqrt(CGAL::squared_distance(points[u], points[v]));
            pq.pop();

            if(distances[local_indices.at(u)][local_indices.at(v)] > t * euc_u_to_v){
                vector<number_t> updatedDistances = Dijkstra_PartitionedCells(points, indices, u, adjMap, local_indices);

                //determine whether the recently calculated paths are less than those stored in distance matrix
                //if so, update distance matrix
                for(size_t i{}; i < n; ++i){

                    number_t newDistance = updatedDistances[i];
                    if(newDistance < distances[local_indices.at(u)][i]){

                        distances[local_indices.at(u)][i] = newDistance;
                        distances[i][local_indices.at(u)] = newDistance;

                    }
                }

                if(updatedDistances[local_indices.at(v)] > t * euc_u_to_v){

                    edges.emplace_back(Edge{u,v});
                    adjMap.at(u).emplace_back(v);
                    adjMap.at(v).emplace_back(u);

                    distances[local_indices.at(u)][local_indices.at(v)] = euc_u_to_v;
                    distances[local_indices.at(v)][local_indices.at(u)] = euc_u_to_v;

                }
            }
        }

//        cout << indices.size() << " | " << edges.size() << endl;

        return CellInfo{distances, local_indices, edges};
    }

    /*
     * Constructs a geometric spanner where all |V| choose 2 pairs (u,v)
     * within the graph contain the invariant that their
     * graphPath(u,v) <= t * d(u,v)
     *
     * Params:
     * points: entire point set
     * t: spanner invariant
     *
     * Returns:
     * Graph consisting of points used, edges created, and adjacency list
     */
    Graph FG_GreedySpanner(const vector<Point>& points,
                           double t){

        size_t n = points.size();

        vector<size_t> indices{};
        indices.reserve(n);
        for(size_t i{}; i < n; ++i){
            indices[i] = i;
        }

        vector<vector<number_t>> distances(n, vector<number_t>{numeric_limits<double>::max()});

        PointPairPQ pairs = sortPairs(points, indices);

        unordered_map<size_t, vector<size_t>> adjMap;
        for(size_t i{}; i < n; ++i){
            adjMap.insert(make_pair(i, vector<size_t>{}));
        }

        vector<Edge> edges{};

        while(!pairs.empty()){
            auto [pts, euc_u_to_v] = pairs.top();
            auto u = pts.first, v = pts.second;
            pairs.pop();

            if(distances[u][v] > euc_u_to_v * t){
                vector<number_t> costs = Dijkstra(points, indices, u, adjMap);

                for(size_t i{}; i < n; ++i){
                    number_t updatedCost = costs[i];
                    if(distances[u][i] > updatedCost){
                        distances[u][i] = updatedCost;
                        distances[i][u] = updatedCost;
                    }
                }

                if(distances[u][v] > euc_u_to_v * t){
                    edges.emplace_back(Edge{u,v});
                    adjMap.at(u).emplace_back(v);
                    adjMap.at(v).emplace_back(u);
                }

            }

        }

        return Graph{points, edges, adjMap};

    }

}

#endif //SP_MER_17_GREEDYSPANNERS_H
