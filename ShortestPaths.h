//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_SHORTESTPATHS_H
#define SP_MER_17_SHORTESTPATHS_H

#include "Utilities.h"
#include "IndexPriorityQueue.h"

namespace spanners{

    /*
     * Used to determine the path costs within a partitioned cell.
     * Specialized to reduce the memory cost as each partition consists of
     * a small subset of points relative to the entire point set
     * As such, a localIndices map is required to translate from
     * original index to local index in order to retrieve values
     *
     * Params:
     * points: entire point set
     * indices: original indices used to index into points
     * u: source vertex
     * adjMap: partitions adjacency list
     * localIndices: map to translate from original indices to local indices
     *
     * Returns:
     * path costs for all vertices from source u
     */
    vector<number_t> Dijkstra_PartitionedCells(const vector<Point> &points,
                                               const vector<size_t> &indices,
                                               size_t u,
                                               const unordered_map<size_t, vector<size_t>> &adjMap,
                                               const unordered_map<size_t, size_t> &local_indices){

        vector<number_t> distances(indices.size(), CGAL_IA_MAX_DOUBLE);
        distances[local_indices.at(u)] = 0;

        IndexPQ::IndexPriorityQueue<size_t, number_t> ipq(indices, distances, IndexPQ::IndexPQType{IndexPQ::IndexPQType::MinHeap});
        ipq.updateKey(u, 0.0);

        while(!ipq.empty()){

            auto [currentNode, cost] = ipq.frontKV();
            ipq.pop();

            for(auto a : adjMap.at(currentNode)){
                if(ipq.contains(a)){

                    number_t distanceFromAdjacent = cost + getDistance(points[currentNode], points[a]);

                    if(distanceFromAdjacent < distances[local_indices.at(a)]){

                        distances[local_indices.at(a)] = distanceFromAdjacent;
                        ipq.updateKey(a, distanceFromAdjacent);

                    }
                }
            }
        }

        return distances;

    }

    /*
     * Traditional Dijkstra
     *
     * Params:
     * points: entire point set
     * indices: used to index points
     * u: source vertex
     * adjMap: graph adjaceny list
     *
     * Returns:
     * path costs for all vertices from source u
     */
    vector<number_t> Dijkstra(const vector<Point> &points,
                              const vector<size_t> &indices,
                              size_t u,
                              const unordered_map<size_t, vector<size_t>> &adjMap){

        vector<number_t> distances(indices.size(), CGAL_IA_MAX_DOUBLE);
        distances[u] = 0;

        IndexPQ::IndexPriorityQueue<size_t, number_t> ipq(indices, distances, IndexPQ::IndexPQType{IndexPQ::IndexPQType::MinHeap});
        ipq.updateKey(u, 0);

        while(!ipq.empty()){

            auto [currentNode, cost] = ipq.frontKV();
            ipq.pop();

            for(auto a : adjMap.at(currentNode)){
                if(ipq.contains(a)){

                    number_t distanceFromAdjacent = cost + getDistance(points[currentNode], points[a]);

                    if(distanceFromAdjacent < distances[a]){

                        distances[a] = distanceFromAdjacent;
                        ipq.updateKey(a, distanceFromAdjacent);

                    }
                }
            }
        }

        return distances;

    }

    /*
     * Traditional a* search algorithm using h = d(i, goal)
     *
     * Params:
     * points: entire point set
     * adjMap: adjacency list
     * u: source
     * v: goal
     *
     * Returns:
     * path cost from source u to goal v
     */
    number_t aStar(const vector<Point> &points,
                   const unordered_map<size_t, vector<size_t>> &adjMap,
                   size_t u, size_t v){

        function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>&)> cmp =
                [](const pair<size_t, number_t> &o1, const pair<size_t, number_t> &o2){
                    return o2.second < o1.second;
                };

        priority_queue<pair<size_t, number_t>, vector<pair<size_t, number_t>>,
                function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>)>> openSet(cmp);

        openSet.push(pair<size_t, number_t>{u, 0});

        unordered_map<size_t, number_t> g_score{}, f_score;
        g_score.insert(pair<size_t, number_t>{u, 0});

        function<number_t(size_t, size_t)> h = [&points](size_t u, size_t v)->number_t {
            return CGAL::sqrt(CGAL::squared_distance(points[u], points[v]));
        };

        f_score.insert(pair<size_t, number_t>{u, h(u,v)});

        unordered_set<size_t> openSetEntries{};
        openSetEntries.insert(u);

        while(!openSet.empty()){

            auto [cur_pt, f] = openSet.top();
            if(cur_pt == v){
                return g_score.at(v);
            }
            openSet.pop();
            openSetEntries.erase(cur_pt);
            if(g_score.find(cur_pt) == g_score.end()){
                g_score.insert(pair<size_t, number_t>{cur_pt, CGAL_IA_MAX_DOUBLE});
            }

            auto adj = adjMap.at(cur_pt);
            for(const auto &a : adj){

                number_t tentative_gScore = g_score.at(cur_pt);

                if(tentative_gScore != CGAL_IA_MAX_DOUBLE){
                    tentative_gScore += CGAL::sqrt(CGAL::squared_distance(points[cur_pt], points[a]));
                }

                if(g_score.find(a) == g_score.end()){
                    g_score.insert(pair<size_t, number_t>{a, CGAL_IA_MAX_DOUBLE});
                }
                if(f_score.find(a) == f_score.end()){
                    f_score.insert(pair<size_t, number_t>{a, CGAL_IA_MAX_DOUBLE});
                }

                if(tentative_gScore < g_score.at(a)){

                    g_score.at(a) = tentative_gScore;

                    f_score.at(a) = tentative_gScore + h(a, v);

                    if(openSetEntries.find(a) == openSetEntries.end()){
                        openSetEntries.insert(a);
                        openSet.push(pair<size_t, number_t>{a, f_score.at(a)});
                    }

                }
            }
        }
        //failed to find goal state
        return -1;
    }
    
}

#endif //SP_MER_17_SHORTESTPATHS_H
