//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_FINDCENTROIDS_H
#define PARTITIONEDSPANNER_FINDCENTROIDS_H

#include "Utilities.h"

namespace spanners{

    vector<Point> findCentroid(const vector<Point> &points,
                               unordered_map<size_t, Aux_QT_Node>& leaf_info,
                               const unordered_map<size_t, Point>& IndexToPoint){
        vector<Point> leaders;

        size_t count{};
        for(auto &i : leaf_info){

            number_t sum_x{}, sum_y{}, avg_x{}, avg_y{};
            for(const auto &pt : i.second.indices){
                sum_x += points[pt].x();
                sum_y += points[pt].y();
            }

            size_t N = i.second.indices.size();
            avg_x = sum_x / N;
            avg_y = sum_y / N;

            Point virtualCentroid{avg_x, avg_y};

            pair<size_t, number_t> closestToVirtual{0, CGAL_IA_MAX_DOUBLE};

            for(const auto &pt : i.second.indices){
                number_t dist = CGAL::sqrt(CGAL::squared_distance(points[pt], virtualCentroid));
                if(dist < closestToVirtual.second){
                    closestToVirtual.first = pt;
                    closestToVirtual.second = dist;
                }
            }

            i.second.leader = closestToVirtual.first;
            i.second.leaderActual = count++;
            leaders.emplace_back(IndexToPoint.at(i.second.leader));

        }

        return leaders;

    }

}

#endif //PARTITIONEDSPANNER_FINDCENTROIDS_H
