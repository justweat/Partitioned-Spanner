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

        auto start = std::chrono::high_resolution_clock::now();

        vector<Point> points = GetPoints();

        double t{};
        cout << "T Val: ";
        cin >> t;

        int cellSize{};
        cout << "Cell Size: ";
        cin >>cellSize;

        Graph spanner = partitionedSpanner(points, cellSize, t);

        auto finish = std::chrono::high_resolution_clock::now();

        cout << "Edges: " << spanner.edges.size() << endl;

        GraphPrinter printer("./", "article");
        printer.autoscale(points.begin(), points.end());
        printer.drawEdges(spanner.edges.begin(), spanner.edges.end(), points);
        printer.display();

    }

}

#endif //SP_MER_17_TESTING_H
