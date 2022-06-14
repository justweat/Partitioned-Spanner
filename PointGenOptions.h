//
// Created by justin on 1/11/22.
//

#ifndef SP_MER_POINTGENOPTIONS_H
#define SP_MER_POINTGENOPTIONS_H

#include "Utilities.h"
#include "PointGenerators.h"

namespace spanners{
    using namespace std;

    vector<Point> GeneratePoints(vector<int> options){

        spanners::RandomPointGenerator_2 pg;
        vector<Point> points{};

        switch(options[0]){
            case 1:
                pg.generatePointsInGalaxy(options[0], options[1], points);
                break;
            case 2:
                pg.generatePointsInsideASquare(options[0], options[1], points);
                break;
            case 3:
                pg.generatePointsInsideASquareNormal(options[0], options[1], points);
                break;
            case 4:
                pg.generatePointsOnASquare(options[0], options[1], points);
                break;
            case 5:
                pg.generateContiguousPointsOnAGrid(options[0], points);
                break;
            case 6:
                pg.generatePointsInsideADisc(options[0], options[1], points);
                break;
            case 7:
                pg.generatePointsOnACircle(options[0], options[1], points);
                break;
            case 8:
                pg.generatePointsOnSpokes(options[0], options[1], points);
                break;
            case 9:
                pg.generateRandomInsideAnnulus(options[0], options[1], options[2], points);
                break;
            case 10:
                pg.generateRandomPointsOnAGrid(options[0], points);
                break;
            default:
                throw logic_error("Invalid point generation option");
        }

        return points;

    }

    vector<Point> GetPoints(){
        vector<Point> points{};
        spanners::RandomPointGenerator_2 pg;
        cout << "1: Galaxy" << endl;
        cout << "2: Inside Square" << endl;
        cout << "3: Inside Square Normal" << endl;
        cout << "4: On Square" << endl;
        cout << "5: Contiguous On Grid" << endl;
        cout << "6: Inside Disc" << endl;
        cout << "7: On Circle" << endl;
        cout << "8: On Spokes" << endl;
        cout << "9: Inside Annulus" << endl;
        cout << "10: On A Grid" << endl;
        int c = 0;
        cout << "Option: ";
        cin >> c;
        int n;
        int k;
        double size, size2;
        switch(c){
            case 1:
                cout << "Num pts: ";
                cin >> n;
                cout << "Num spokes: ";
                cin >> k;
                pg.generatePointsInGalaxy(n, k, points);
                break;
            case 2:
                cout << "Num pts: ";
                cin >> n;
                cout << "Size of square: ";
                cin >> size;
                pg.generatePointsInsideASquare(n, size, points);
                break;
            case 3:
                cout << "Num of pts per cluster: ";
                cin >> n;
                cout << "Num of clusters: ";
                cin >> k;
                pg.generatePointsInsideASquareNormal(n, k, points);
                break;
            case 4:
                cout << "Num pts: ";
                cin >> n;
                cout << "Size of square: ";
                cin >> size;
                pg.generatePointsOnASquare(n, size, points);
                break;
            case 5:
                cout << "Num pts: ";
                cin >> n;
                pg.generateContiguousPointsOnAGrid(n, points);
                break;
            case 6:
                cout << "Num pts: ";
                cin >> n;
                cout << "Radius: ";
                cin >> size;
                pg.generatePointsInsideADisc(n, size, points);
                break;
            case 7:
                cout << "Num pts: ";
                cin >> n;
                cout << "Size of square: ";
                cin >> size;
                pg.generatePointsOnACircle(n, size, points);
                break;
            case 8:
                cout << "Num pts: ";
                cin >> n;
                cout << "Num spokes: ";
                cin >> k;
                pg.generatePointsOnSpokes(n, k, points);
                break;
            case 9:
                cout << "Num pts: ";
                cin >> n;
                cout << "Radius 1: ";
                cin >> size;
                cout << "Radius 2: ";
                cin >> size2;
                pg.generateRandomInsideAnnulus(n, size, size2, points);
                break;
            case 10:
                cout << "Num pts: ";
                cin >> n;
                pg.generateRandomPointsOnAGrid(n, points);
                break;
            default:
                throw logic_error("Invalid point generation option");
        }
        return points;
    }
}

#endif //SP_MER_POINTGENOPTIONS_H
