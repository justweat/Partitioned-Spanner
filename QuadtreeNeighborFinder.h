//
// Created by justin on 2/27/22.
//

#ifndef SP_MER_17_QUADTREENEIGHBORFINDER_H
#define SP_MER_17_QUADTREENEIGHBORFINDER_H

#include "Utilities.h"

namespace spanners{

    /*
     * Recurses quadtree nodes until reaching leaf
     * Used for: S node finding N neighbors
     */
    void findLeafNeighbors_North(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_North(cur[0], neighbors);
        findLeafNeighbors_North(cur[1], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: E node finding W neighbors
    */
    void findLeafNeighbors_West(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_West(cur[3], neighbors);
        findLeafNeighbors_West(cur[1], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: N node finding S neighbors
    */
    void findLeafNeighbors_South(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_South(cur[2], neighbors);
        findLeafNeighbors_South(cur[3], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: NW node finding SE neighbors
    */
    void findLeafNeighbors_SE(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_SE(cur[2], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: NE node finding SW neighbors
    */
    void findLeafNeighbors_SW(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_SW(cur[3], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: W node finding E neighbors
    */
    void findLeafNeighbors_East(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_East(cur[0], neighbors);
        findLeafNeighbors_East(cur[2], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: SW node finding NE neighbors
    */
    void findLeafNeighbors_NE(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_NE(cur[0], neighbors);
    }

    /*
    * Recurses quadtree nodes until reaching leaf
    * Used for: SE node finding NW neighbors
    */
    void findLeafNeighbors_NW(QT_Node cur, vector<QT_Node> &neighbors){
        if(cur.is_leaf()){
            neighbors.emplace_back(cur);
            return;
        }
        findLeafNeighbors_NW(cur[1], neighbors);
    }

}



#endif //SP_MER_17_QUADTREENEIGHBORFINDER_H
