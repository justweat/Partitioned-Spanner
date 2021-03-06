//
// Created by justin on 5/16/22.
//

#include "QuadtreeNeighborFinder.h"
#include "ResolveContiguousNeighbors.h"
#include "ResolveDistantPairs.h"
#include "InitializePartitions.h"

#ifndef SP_MER_17_PARTITIONEDSPANNER_H
#define SP_MER_17_PARTITIONEDSPANNER_H

namespace spanners{

    /*
     * Geometric spanner
     * -----------------
     * Defined as a graph with an invariant t >= 1 where
     * all |V| choose 2 pairs of points consist of paths
     * such that path(u, v) <= d(u, v) * t
     */

    /*
     * Bridging heuristic
     * ------------------
     * Contiguous neighbors utilize the following approach:
     * After the partitions have been initialized,
     * the current state of the graph consists of disjoint cells.
     * To merge contiguous neighbors (Partition_i, Partition_j),
     * sort pairs of points between each partition from i to j in non-decreasing order denoted (Pi, Pj)
     * based on their Euclidean distance.
     * For each bridge in existence where Bi is the entrance from Partition_i and Bj is the exit into Partition_j,
     * reuse the distance matrices from the partition initialization and
     * determine whether there exists a viable path such that
     * Path(Pi, Bi) + d(Bi, Bj) + Path(Bj, Pj) <= t * d(Pi, Pj)
     * If true, move to next pair.
     * If false, create new edge and add this to Bridges.
     *
     * Distant pair resolution consists of the same approach except we add two initial bridges:
     * The path from the closest pairs (Pi, Pj) that traverses the spanner in its current state
     * The leader spanner path from the centroids of Partition_i to Partition_j
     */

    //TODO: find more elegant way to structure the algo without the need of a default boolean
    /*
    * Constructs a geometric spanner
    * Accomplished by partitioning the points into [1, cellSize] disjoint cells,
    * then merging those cells to maintain the invariant using the Bridging heuristic described above
    * all while taking advantage of the concurrent nature of this D&C algorithm
    * and reusing initial calculations for subsequent resolution stages
    *
    * Params:
    * points: G.V
    * cellSize: partition maximum
    * leaderSpannerConstructor: initial call for constructor should be false by default
    * numOfThreads: max number of threads used
    *
    * Returns:
    * Graph consisting of points, edges, and adjacency list
    */
    Graph partitionedSpanner(vector<Point> &points,
                             size_t cellSize,
                             double t,
                             unsigned int numOfThreads = 1){

//        cout << "Spanner for " << points.size() << endl;

        size_t n = points.size();

        /*
         * Base case for recursive leaderSpanner calls
         */
        if(n <= cellSize){
            return FG_GreedySpanner(points, t);
        }

        numOfThreads = min(numOfThreads, thread::hardware_concurrency());

        if(numOfThreads < 1) numOfThreads = 1;

        Quadtree qt(points);
        qt.refine(numeric_limits<size_t>::max(), cellSize);

        unordered_map<Point, size_t> PointToIndex{};
        unordered_map<size_t, Point> IndexToPoint{};

        for(size_t i{}; i < n; ++i){
            PointToIndex.insert(pair<Point, size_t>{points[i], i});
            IndexToPoint.insert(pair<size_t, Point>{i, points[i]});
        }

        //Unique ptrs to all partitions
        vector<unique_ptr<Partition>> leaves{};

        //Used for finding continguous neighbors
        map<QT_Node, size_t> leafIdentifier{};

        /*
         * Decompose the point set into partitions while
         * initializing the vars required to accomplish the merging procedures
         */
        size_t partitionIndex{};
        for(Quadtree::Node i : qt.traverse<CGAL::Orthtrees::Leaves_traversal>()){
            if(i.empty()) continue;
            vector<size_t> indices;
            vector<Point> pts;
            for(const auto &p : i){
                indices.emplace_back(PointToIndex.at(p));
                pts.emplace_back(p);
            }
            Partition auxNode{indices, partitionIndex, pts};
            leaves.push_back(make_unique<Partition>(auxNode));
            leafIdentifier.insert(pair<QT_Node, size_t>{i, partitionIndex});
            ++partitionIndex;
        }

        /*
         * Launch thread(s) responsible for establishing the initial state of each partition
         * For each partition:
         * Construct greedy spanner within the cell
         * Find bounding box
         * Find centroid
         */
        thread iniPar{initializePartitions, ref(leaves), points, t, numOfThreads};

        /*
         * For each partition, search for all contiguous neighbors
         * Contiguous neighbor in this scenario defined as all
         * North, East, NorthEast, and NorthWest partitions
         * adjacent to the ith partition
         */
        ContiguousNeighbors contiguousNeighbors = findContiguousNeighbors(leafIdentifier);

        /*
         * Join partition initialization here
         */
        iniPar.join();

        //Safety check to ensure that cell partitions correctly initialized
//        for(const auto& leaf : leaves){
//            for(const auto& i : leaf->distances){
//                for(const auto& j : i){
//                    assert(j != CGAL_IA_MAX_DOUBLE && "Partition greedy construction failed");
//                }
//            }
//        }

        /*
         * Launch thread to merge contiguous partition neighbors as defined above
         */
        thread resConNei{resolveContiguousNeighbors, cref(points), t, ref(leaves), cref(contiguousNeighbors), numOfThreads};

        /*
         * Use centroids to establish a "leader spanner"
         * Required to ensure that partitions are connected
         */

        vector<Point> leaderPoints{};
        //used to map leader spanner edges back to original spanner and update adjacency list
        unordered_map<size_t, size_t> leaderPointTranslator{};
        size_t leaderIndex{};

        for(const auto& leaf : leaves){
            size_t originalIndex = leaf->leader;
            size_t newIndex = leaderIndex++;
            leaderPointTranslator.insert(pair<size_t, size_t>{newIndex, originalIndex});
            leaderPoints.push_back(points[leaf->leader]);
//            leaderPoints.push_back(IndexToPoint.at(leaf->leader));
        }

        /*
         * Construct leader spanner
         */
        Graph leaderSpanner = partitionedSpanner(leaderPoints, cellSize, t, numOfThreads);
//        Graph leaderSpanner = FG_GreedySpanner(leaderPoints, t);

        /*
         * Initialize adjacency map
         */
        unordered_map<size_t, vector<size_t>> adjMap;
        for(size_t i{}; i < n; ++i){
            adjMap.insert(make_pair(i, vector<size_t>{}));
        }

        /*
         * Join contiguous neighbor resolution thread here
         */
        resConNei.join();

        /*
         * Ensure adjacency map reflects the current state of the spanner
         * after constructing partitions and merging their contiguous neighbors
         * Required for distant pair resolution described below
         */
        for(const auto& leaf : leaves){

            for(const auto& e : leaf->greedyEdges){
                adjMap.at(e.first).push_back(e.second);
                adjMap.at(e.second).push_back(e.first);
            }

            for(const auto& e : leaf->contiguousResolutionEdges){
                adjMap.at(e.first).push_back(e.second);
                adjMap.at(e.second).push_back(e.first);
            }

        }

        for(const auto& e : leaderSpanner.edges){
            adjMap.at(leaderPointTranslator.at(e.first)).push_back(leaderPointTranslator.at(e.second));
            adjMap.at(leaderPointTranslator.at(e.second)).push_back(leaderPointTranslator.at(e.first));
        }

        /*
         * For each partition,
         * find all candidate distant pairs such that they were
         * not used in contiguous neighbor resolution and
         * satisfy the equation described at the function implementation
         * For each candidate pair,
         * concurrently merge each partition pair using the distant pair Bridging heuristic
         */
        vector<Edge> distantEdges = resolveDistantPairs(points,
                                                        adjMap,
                                                        leaderPoints,
                                                        leaves,
                                                        t,
                                                        leaderSpanner,
                                                        contiguousNeighbors.contiguousNeighborCheck);

        /*
         * Spanner construction complete
         * Merge all edges created throughout and return final graph
         */
        vector<Edge> totalEdges{};

        size_t sumGre{}, sumRes{};

        for(const auto& leaf : leaves){

            for(const auto& e : leaf->greedyEdges){
                ++sumGre;
                totalEdges.push_back(e);
            }

            for(const auto& e : leaf->contiguousResolutionEdges){
                ++sumRes;
                totalEdges.push_back(e);
            }

        }

        for(const auto& e : distantEdges){
            totalEdges.push_back(e);
        }

        for(auto& e : leaderSpanner.edges){
            totalEdges.emplace_back(leaderPointTranslator.at(e.first), leaderPointTranslator.at(e.second));
        }

        cout << "Greedy Edge Count: " << sumGre << endl;
        cout << "Contiguous Resolution Edge Count: " << sumRes << endl;
        cout << "Distant Resolution Edge Count: " << distantEdges.size() << endl;
        cout << "Leader Spanner Edge Count: " << leaderSpanner.edges.size() << endl;

        return Graph{points, totalEdges, adjMap};

    }

}

#endif //SP_MER_17_PARTITIONEDSPANNER_H
