//
// Created by justin on 6/21/22.
//

#ifndef FAST_SPARSE_SPANNER_INITIALIZEPARTITIONS_H
#define FAST_SPARSE_SPANNER_INITIALIZEPARTITIONS_H

#include "GreedySpanners.h"
#include "FindCentroids.h"

#include <thread>

namespace spanners{

    /*
     * Initializes the partitions to consist of:
     * -Greedy spanner
     * -Bounding box
     * -Leader/centroid point
     *
     * Params:
     * points: entire point set
     * partitions: unique pts to each partition
     * t: spanner invariant
     * begin/end: indices within partitions for this thread to initialize [begin, end]
     *
     * Returns:
     * void (results stored using partition ptr)
     */
    void createPartitions(const vector<Point>& points,
                          vector<unique_ptr<Partition>>& partitions,
                          double t,
                          size_t begin,
                          size_t end){

        for(size_t i {begin}; i < end; ++i){

            unique_ptr<Partition>& currentPartition = partitions[i];

            CellInfo cell = FG_GreedySpanner_PartitionedCell(points, currentPartition->indices, currentPartition->adjMap, t);

            currentPartition->greedyEdges = move(cell.partitionEdges);
            currentPartition->distances = move(cell.partitionDistanceMatrix);
            currentPartition->localIndices = move(cell.partitionLocalIndices);

            currentPartition->bounding_box = CGAL::bbox_2(currentPartition->points.begin(), currentPartition->points.end());

            currentPartition->leader = findCentroid(points, currentPartition);

        }

    }

    /*
     * Launch threads to initialize partitions
     *
     * partitions: unique ptrs to each leaf/partition
     * points: entire point set
     * t: spanner invariant
     * numOfThreads: max number of threads
     *
     * Returns:
     * void (results stored using partition ptr)
     */
    void initializePartitions(vector<unique_ptr<Partition>>& partitions,
                                      const vector<Point>& points,
                                      double t,
                                      unsigned int numOfThreads = 1){

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

        size_t n = partitions.size();

        size_t block = n / numOfThreads;
        if(block == 0){
            block = n;
        }

        vector<thread> threads{};

        for(size_t i{}; i < n; i += block){
            threads.emplace_back(createPartitions, cref(points), ref(partitions), t, i, min(n - 1, i + block));
        }

        for(auto& thread : threads){
            thread.join();
        }

    }

}

#endif //FAST_SPARSE_SPANNER_INITIALIZEPARTITIONS_H
