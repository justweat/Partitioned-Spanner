//
// Created by justin on 6/21/22.
//

#ifndef FAST_SPARSE_SPANNER_STRETCHFACTOR_H
#define FAST_SPARSE_SPANNER_STRETCHFACTOR_H

#include "Utilities.h"
#include "ShortestPaths.h"

#include <thread>

namespace spanners{

    struct StretchFactorResult{
        vector<Edge> errors{};
        number_t stretchFactor{};
    };

    void MeasureStretchFactor(const vector<Point>& points,
                              const vector<size_t>& indices,
                              const unordered_map<size_t, vector<size_t>>& adjMap,
                              number_t t,
                              size_t begin,
                              size_t end,
                              StretchFactorResult& results,
                              mutex& resultsLock){

        vector<Edge> errorPaths{};
        number_t threadStretchFactor = numeric_limits<double>::lowest();

        size_t n = points.size();

        for(size_t i{begin}; i <= end; ++i){

            vector<number_t> costs = Dijkstra(points, indices, i, adjMap);

            for(size_t j{}; j < i; ++j){
                costs[j] /= getDistance(points[i], points[j]);
            }

            for(size_t j{i + 1}; j < n; ++j){
                costs[j] /= getDistance(points[i], points[j]);
            }

            for(size_t j{}; j < n; ++j){

                threadStretchFactor = max(threadStretchFactor, costs[j]);

                if(costs[j] > t){
                    errorPaths.emplace_back(Edge{i, j});
                }

            }
        }

        lock_guard<mutex> lk(resultsLock);
        results.stretchFactor = max(results.stretchFactor, threadStretchFactor);
        for(auto e : errorPaths){
            results.errors.emplace_back(move(e));
        }

    }

    StretchFactorResult ParallelStretchFactor(const vector<Point>& points,
                                              const unordered_map<size_t, vector<size_t>>& adjMap,
                                              number_t t,
                                              unsigned int numOfThreads = 1){

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

        size_t n = points.size();

        size_t block = n / numOfThreads;
        if(block == 0){
            block = 1;
        }

        StretchFactorResult results{};
        mutex resultsLock{};

        vector<size_t> indices{};
        for(size_t i{}; i < n; ++i){
            indices.emplace_back(i);
        }

        vector<thread> threads{};
        for(size_t i{}; i < n; i += block){
            threads.emplace_back(MeasureStretchFactor, cref(points), cref(indices), cref(adjMap), t, i, min(n - 1, i + block), ref(results), ref(resultsLock));
        }

        for(auto& thread : threads){
            thread.join();
        }

        return results;

    }

}

#endif //FAST_SPARSE_SPANNER_STRETCHFACTOR_H
