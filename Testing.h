//
// Created by justin on 2/28/22.
//

#ifndef SP_MER_17_TESTING_H
#define SP_MER_17_TESTING_H

#include <CGAL/Quadtree.h>
#include <CGAL/Orthtree/Traversals.h>
#include <chrono>
#include <utility>

#include "PointGenOptions.h"
#include "GraphPrinter.h"
#include "PartitionedSpanner.h"
#include "StretchFactor.h"



namespace spanners{

    void testingPartitionedSpannerAutomated(double t, int cellSize, bool print, vector<int> options){

        auto start = std::chrono::high_resolution_clock::now();

        vector<Point> points = GeneratePoints(std::move(options));

        Graph spanner = partitionedSpanner(points, cellSize, t);

        auto finish = std::chrono::high_resolution_clock::now();

        cout << "Edges: " << spanner.edges.size() << endl;

        if(print){
            GraphPrinter printer("./", "article");
            printer.autoscale(points.begin(), points.end());
            printer.drawEdges(spanner.edges.begin(), spanner.edges.end(), points);
            printer.display();
        }

    }

    void testingPartitionSpanner(){

        vector<Point> points = GetPoints();

        double t{};
        cout << "T Val: ";
        cin >> t;

        int cellSize{};
        cout << "Cell Size: ";
        cin >>cellSize;

        unsigned int numOfThreads{};
        cout << "Num of threads: ";
        cin >> numOfThreads;

        bool sf{};
        cout << "Measure SF: ";
        cin >> sf;

        bool print{};
        cout <<"Print: ";
        cin >> print;

        auto start = std::chrono::high_resolution_clock::now();

        Graph spanner = partitionedSpanner(points, cellSize, t, false, numOfThreads);

        auto finish = std::chrono::high_resolution_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::microseconds>(finish-start);

        cout << "Time: " << (time.count()/1000000)/60 << "m " << (time.count()/1000000)%60 << "s" << endl;
        cout << "Edges: " << spanner.edges.size() << endl;
        cout << "Points: " << spanner.points.size() << endl;
        cout << "M/N: " << static_cast<double>(spanner.edges.size()) / spanner.points.size() << endl;

        StretchFactorResult results;
        if(sf){
            results = ParallelStretchFactor(points, spanner.adjacencyMap, t, numOfThreads);
            cout << "SF: " << results.stretchFactor << "\n";
            cout << "Red Edges: " << results.errors.size() / 2 << "\n";
        }

        if(print && spanner.edges.size() + results.errors.size() < 50000){
            GraphPrinter printer("./", "article");
            printer.autoscale(points.begin(), points.end());
            printer.drawEdges(spanner.edges.begin(), spanner.edges.end(), points);
            printer.drawEdges(results.errors.begin(), results.errors.end(), points, printer.redEdge);
            printer.display();
        }

    }

}

#endif //SP_MER_17_TESTING_H
