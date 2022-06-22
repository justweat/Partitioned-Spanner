//
// Created by justin on 6/14/22.
//

#ifndef PARTITIONEDSPANNER_FINDCENTROIDS_H
#define PARTITIONEDSPANNER_FINDCENTROIDS_H

#include "Utilities.h"

namespace spanners{

    size_t findCentroid(const vector<Point> &points,
                        unique_ptr<Partition>& partition){

        number_t sum_x{}, sum_y{};
        for(const auto &pt : partition->indices){
            sum_x += points[pt].x();
            sum_y += points[pt].y();
        }

        size_t n = partition->indices.size();
        number_t avg_x = sum_x / n;
        number_t avg_y = sum_y / n;

        Point virtualCentroid{avg_x, avg_y};

        pair<size_t, number_t> closestToVirtual{0, CGAL_IA_MAX_DOUBLE};

        for(const auto &pt : partition->indices){
            number_t dist = CGAL::sqrt(CGAL::squared_distance(points[pt], virtualCentroid));
            if(dist < closestToVirtual.second){
                closestToVirtual.first = pt;
                closestToVirtual.second = dist;
            }
        }

        return closestToVirtual.first;

    }

}

#endif //PARTITIONEDSPANNER_FINDCENTROIDS_H
