//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_SHORTESTPATHS_H
#define SP_MER_17_SHORTESTPATHS_H

#include "Utilities.h"
#include "IndexPriorityQueue.h"

namespace spanners{

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

    number_t DijkstaReturnDistances_PartitionedCells(
                                    const vector<Point> &points,
                                    const vector<size_t> &indices,
                                    size_t u, size_t v,
                                    const unordered_map<size_t, vector<size_t>> &adjMap,
                                    const unordered_map<size_t, size_t> &local_indices,
                                    vector<number_t> &distances){

        size_t n = indices.size();

        vector<number_t> local_distances(n, CGAL_IA_MAX_DOUBLE);
        local_distances[local_indices.at(u)] = 0;


        function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>&)> cmp =
                [](const pair<size_t, number_t> &o1, const pair<size_t, number_t> &o2){
                    return o2.second < o1.second;
                };

        priority_queue<pair<size_t, number_t>, vector<pair<size_t, number_t>>,
                function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>)>> pq(cmp);

        pq.push(pair<size_t, number_t>{u, 0});

        while(!pq.empty()){

            auto [cur_point, path_distance] = pq.top();
            pq.pop();

            auto adj_nodes = adjMap.at(cur_point);

            for(const auto &i : adj_nodes){
                number_t dist_curToAdj = local_distances[local_indices.at(cur_point)] + CGAL::sqrt(CGAL::squared_distance(points[cur_point], points[i]));
                if(local_distances[local_indices.at(i)] > dist_curToAdj && local_distances[local_indices.at(cur_point)] != CGAL_IA_MAX_DOUBLE){
                    local_distances[local_indices.at(i)] = dist_curToAdj;
                    pq.push(pair<size_t, number_t>{i, dist_curToAdj});
                }
            }
        }

        for(const auto &i : indices){
            distances[local_indices.at(i)] = local_distances[local_indices.at(i)];
        }

        return local_distances[local_indices.at(v)];

    }

    number_t Dijksta_RestrictedPath_ReturnDistanceVector(const vector<Point> &points,
                                    const vector<size_t> &indices,
                                    size_t u, size_t v,
                                    const unordered_map<size_t, vector<size_t>> &adjMap,
                                    const unordered_map<size_t, size_t> &local_indices,
                                    const unordered_set<size_t> &permitted_indices,
                                    vector<number_t> &distances){

        size_t n = indices.size();

        vector<number_t> local_distances(n, CGAL_IA_MAX_DOUBLE);
        local_distances[local_indices.at(u)] = 0;


        function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>&)> cmp =
                [](const pair<size_t, number_t> &o1, const pair<size_t, number_t> &o2){
                    return o2.second < o1.second;
                };

        priority_queue<pair<size_t, number_t>, vector<pair<size_t, number_t>>,
                function<bool(const pair<size_t, number_t>&, const pair<size_t, number_t>)>> pq(cmp);

        pq.push(pair<size_t, number_t>{u, 0});

        while(!pq.empty()){

            auto [cur_point, path_distance] = pq.top();
            pq.pop();

            auto adj_nodes = adjMap.at(cur_point);

            for(const auto &i : adj_nodes){
                if(permitted_indices.find(i) != permitted_indices.end()){
                    number_t dist_curToAdj = path_distance + CGAL::sqrt(CGAL::squared_distance(points[cur_point], points[i]));
                    if(local_distances[local_indices.at(i)] > dist_curToAdj){
                        local_distances[local_indices.at(i)] = dist_curToAdj;
                        pq.push(pair<size_t, number_t>{i, dist_curToAdj});
                    }
                }
            }
        }

        for(const auto &i : indices){
            distances[local_indices.at(i)] = local_distances[local_indices.at(i)];
        }

        return local_distances[local_indices.at(v)];

    }

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

//        unordered_map<size_t, size_t> cameFrom{};

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
//                number_t path{};
//                size_t traveler = v;
//                while(traveler != u){
//                    path += CGAL::sqrt(CGAL::squared_distance(points[traveler], points[cameFrom.at(traveler)]));
//                    traveler = cameFrom.at(traveler);
//                }
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

//                    if(cameFrom.find(a) == cameFrom.end()){
//                        cameFrom.insert(pair<size_t, size_t>{a, cur_pt});
//                    }else{
//                        cameFrom.at(a) = cur_pt;
//                    }

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
