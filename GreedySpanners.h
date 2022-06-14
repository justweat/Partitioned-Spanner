//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_GREEDYSPANNERS_H
#define SP_MER_17_GREEDYSPANNERS_H

#include "Utilities.h"
#include "ShortestPaths.h"
#include "SortPairs.h"
#include "CGAL/Real_timer.h"


namespace spanners{

    CellInfo FG_GreedySpannerPartitionedCell(
                                const vector<Point> &points,
                                const vector<size_t> &indices,
                                unordered_map<size_t, vector<size_t>> &adjMap,
                                double t
                                ){

        vector<Edge> edges{};
        size_t n = indices.size();

        unordered_map<size_t, size_t> local_indices{};
        for(size_t i{}; i < n; ++i){
            local_indices.insert(make_pair(indices[i], i));
        }

        PointPairPQ pq = sortPairs(points, indices);

        vector<vector<number_t>> distances(n, vector<number_t>(n, CGAL_IA_MAX_DOUBLE));

        while(!pq.empty()){

            auto [pts, euc_u_to_v] = pq.top();
            auto u = pts.first, v = pts.second;
            euc_u_to_v = CGAL::sqrt(CGAL::squared_distance(points[u], points[v]));
            pq.pop();

            if(distances[local_indices.at(u)][local_indices.at(v)] > t * euc_u_to_v){
                number_t updated_dist = DijkstaReturnDistances_PartitionedCells(points, indices, u, v, adjMap, local_indices, distances[local_indices.at(u)]);
                //updated both distance vectors to min for dist(u,v) and dist(v,u)
                if(updated_dist > t * euc_u_to_v){
                    edges.emplace_back(Edge{u,v});
                    adjMap.at(u).emplace_back(v);
                    adjMap.at(v).emplace_back(u);
                    distances[local_indices.at(u)][local_indices.at(v)] = euc_u_to_v;
                    distances[local_indices.at(v)][local_indices.at(u)] = euc_u_to_v;
                }
            }
        }
        return CellInfo{distances, local_indices, edges};
    }

}

#endif //SP_MER_17_GREEDYSPANNERS_H
