#include "Testing.h"

using namespace std;

int main(int argc, char** argv) {

    //argv:
    //print, tVal, cellSize, point generation options

    //point generation options
    //numPts, spannerType, spannerTypeOptions

    if(argc == 1){
        spanners::testingPartitionSpanner();
    }else{
        bool print = stoi(argv[1]);

        double t = stod(argv[2]);

        int cellSize = stoi(argv[3]);

        vector<int> options{};
        for(int i = 3; i < argc; ++i){
            options.emplace_back(stoi(argv[i]));
        }

        switch(options[0]){
            case 5:
            case 10:
                assert(options.size() == 1);
                break;
            case 1:
            case 2:
            case 3:
            case 4:
            case 6:
            case 7:
            case 8:
                assert(options.size() == 2);
                break;
            case 9:
                assert(options.size() == 3);
                break;
            default:
                throw logic_error("Invalid point generation option");
        }

        spanners::testingPartitionedSpannerAutomated(t, cellSize, print, options);
    }

    return 0;
}
