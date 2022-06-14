//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H
#define PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H

#include "Utilities.h"
#include "QuadtreeNeighborFinder.h"
#include "SortPairs.h"

namespace spanners{

    vector<Edge> resolveContiguousNeighbors(const vector<Point> &points,
                                            double t,
                                            QT_Node node,
                                            unordered_map<size_t, Aux_QT_Node>& leaf_info,
                                            const map<QT_Node, size_t>& leaf_identifier,
                                            Direction direction,
                                            set<pair<size_t, size_t>>& contig_neighbor_check){

        vector<Edge> resolved_edges{};

        QT_Node adjacent;
        if(direction == Direction::North){

            adjacent = node.adjacent_node(11);

        }else if(direction == Direction::East){

            adjacent = node.adjacent_node(01);

        }else if(direction == Direction::NorthEast){

            adjacent = node.adjacent_node(01);
            if(!adjacent.is_null()){
                adjacent = adjacent.adjacent_node(11);
            }

        }else if(direction == Direction::NorthWest){

            adjacent = node.adjacent_node(00);
            if(!adjacent.is_null()){
                adjacent = adjacent.adjacent_node(11);
            }

        }else{
            throw logic_error("Not a valid direction for contiguous neighbor resolution");
        }

        if(!adjacent.is_null()){

            auto* aux = &leaf_info.at(leaf_identifier.at(node));
            vector<QT_Node> neighbors;
            if(direction == Direction::North){

                findLeafNeighbors_North(adjacent, neighbors);

            }else if(direction == Direction::East){

                findLeafNeighbors_East(adjacent, neighbors);

            }else if(direction == Direction::NorthEast){

                findLeafNeighbors_NE(adjacent, neighbors);

            }else{

                findLeafNeighbors_NW(adjacent, neighbors);

            }

            vector<size_t> indices_of_neighbors{};
            unordered_map<size_t, size_t> neighbor_map{};
            vector<pair<vector<vector<number_t>>*, unordered_map<size_t, size_t>*>> bridges;
            vector<pair<size_t, size_t>> connections{};
            size_t index{};
            for(const auto &j : neighbors){
                if(j.size() == 0) continue;
                auto* bot = &leaf_info.at(leaf_identifier.at(j));
                for(const auto &k : bot->indices){
                    indices_of_neighbors.emplace_back(k);
                    neighbor_map.insert(pair<size_t, size_t>{k, index});
                }
                contig_neighbor_check.insert(pair<size_t, size_t>{min(aux->index, bot->index), max(aux->index, bot->index)});
                bridges.emplace_back(make_pair(&bot->distances, &bot->local_indices));
                ++index;
            }
            PointPairPQ pairs = sortPairs_DisjointCells(points, aux->indices, indices_of_neighbors);
            while(!pairs.empty()){
                auto [pts, dist] = pairs.top();
                pairs.pop();
                size_t u = pts.first, v = pts.second;

                bool found_path{};
                size_t in = neighbor_map.at(v);
                for(const auto &c : connections){
                    if(neighbor_map.at(c.second) == in){

                        vector<vector<number_t>>& ref = *bridges[in].first;
                        unordered_map<size_t, size_t>& ind = *bridges[in].second;

                        if(ref[ind.at(c.second)][ind.at(v)] +
                           aux->distances[aux->local_indices.at(c.first)][aux->local_indices.at(u)] +
                           CGAL::sqrt(CGAL::squared_distance(points[c.first], points[c.second])) <
                           t * CGAL::sqrt(CGAL::squared_distance(points[u], points[v]))){
                            found_path = true;
                            break;
                        }
                    }
                }

                if(!found_path){
                    connections.emplace_back(pair<size_t,size_t>{u,v});
                    resolved_edges.emplace_back(Edge{u,v});
                }
            }
        }

        return resolved_edges;

    }

}

#endif //PARTITIONEDSPANNER_RESOLVECONTIGUOUSNEIGHBORS_H
