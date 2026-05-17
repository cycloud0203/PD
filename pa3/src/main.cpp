#define _GLIBCXX_USE_CXX11_ABI 0
#include "Placement.h"
#include "Util.h"
#include "GlobalPlacer.h"
#include "arghandler.h"
#include "DPlace.h"
#include "TetrisLegal.h"
#include "ParamPlacement.h"

#include <cstring>
#include <iostream>
#include <string>
#include <time.h>

#include <omp.h>

using namespace std;

bool handleArgument(const int& argc, char* argv[], CParamPlacement& param) {
    int argumentIndex = 0;
    if (argc > 2 && strcmp(argv[1] + 1, "aux") == 0) {
        param.auxFilename = string(argv[2]);
        argumentIndex = 3;
    } else {
        cout << "Usage: " << argv[0] << " -aux benchmark.aux" << endl;
        return false;
    }

    while (argumentIndex < argc) {
        if (strlen(argv[argumentIndex]) <= 1) {
            ++argumentIndex;
            continue;
        }

        if (strcmp(argv[argumentIndex] + 1, "noglobal") == 0) {
            param.bRunGlobal = false;
        } else if (strcmp(argv[argumentIndex] + 1, "nolegal") == 0) {
            param.bRunLegal = false;
        } else if (strcmp(argv[argumentIndex] + 1, "nodetail") == 0) {
            param.bRunDetail = false;
        } else if (strcmp(argv[argumentIndex] + 1, "loadpl") == 0) {
            param.plFilename = string(argv[++argumentIndex]);
        }
        ++argumentIndex;
    }
    return true;
}

int main(int argc, char* argv[]) {
    gArg.Init(argc, argv);
    if (!handleArgument(argc, argv, param)) {
        return -1;
    }

    Placement placement;
    placement.readBookshelfFormat(param.auxFilename, param.plFilename);

    cout << "Benchmark: " << placement.name() << endl;
    cout << format("HPWL: %.f", placement.computeHpwl()) << endl;
    cout << format("Memory usage: %.1f MB", getCurrentMemoryUsage()) << endl;
    cout << format("Core region: (%.f,%.f)-(%.f,%.f)",
                    placement.boundryLeft(),
                    placement.boundryBottom(),
                    placement.boundryRight(),
                    placement.boundryTop())
         << endl;

    double orig_wirelength = 0.0;
    double gp_wirelength = 0.0;
    double lg_wirelength = 0.0;
    double dp_wirelength = 0.0;
    bool bLegal = false;

    time_t total_time = 0;
    time_t global_time_start = time(NULL);
    time_t total_global_time = 0;

    if (param.bRunGlobal) {
        cout << endl << "////// Global Placement ///////" << endl;

        constexpr int kNumThreads = 4;
        omp_set_dynamic(0);
        omp_set_num_threads(kNumThreads);
        omp_set_schedule(omp_sched_static, 0);
        cout << "Using " << kNumThreads << " OpenMP threads" << endl;

        GlobalPlacer globalPlacer;
        globalPlacer.setup(placement);
        globalPlacer.configureBenchmark();

        for (int argumentIndex = 1; argumentIndex < argc; ++argumentIndex) {
            if (strcmp(argv[argumentIndex], "-seed") == 0 && argumentIndex + 1 < argc) {
                globalPlacer.seed = std::atoi(argv[argumentIndex + 1]);
                break;
            }
        }

        cout << "Set seed: " << globalPlacer.seed << endl;
        globalPlacer.solve(globalPlacer.seed);
        globalPlacer.exportDetailedPlot("init.plt", false);

        placement.outputBookshelfFormat(placement.name() + ".gp.pl");

        gp_wirelength = placement.computeHpwl();
        printf("\nHPWL: %.0f\n", gp_wirelength);
        total_global_time = time(NULL) - global_time_start;
        total_time += total_global_time;
    }

    time_t legal_time_start = time(NULL);
    time_t total_legal_time = 0;
    if (param.bRunLegal) {
        cout << endl << "////// Legalization ///////" << endl;
        orig_wirelength = placement.computeHpwl();

        CTetrisLegal legal(placement);
        bLegal = legal.Solve(0.8);
        if (bLegal) {
            cout << "legalization success!" << endl;
        } else {
            cout << "legalization fail!" << endl;
        }

        placement.outputBookshelfFormat(placement.name() + ".lg.pl");

        lg_wirelength = placement.computeHpwl();
        printf("\nHPWL: %.0f (%3.2f%%)\n",
               lg_wirelength,
               ((lg_wirelength - orig_wirelength) / orig_wirelength) * 100.0);
        total_legal_time = time(NULL) - legal_time_start;
        total_time += total_legal_time;
    }

    time_t detail_time_start = time(NULL);
    time_t total_detail_time = 0;
    if (param.bRunDetail && bLegal) {
        cout << endl << "////// Detail Placement ///////" << endl;
        orig_wirelength = placement.computeHpwl();

        CDetailPlacer dplacer(placement);
        dplacer.DetailPlace();

        placement.outputBookshelfFormat(placement.name() + ".dp.pl");

        dp_wirelength = placement.computeHpwl();
        printf("\nHPWL: %.0f (%3.2f%%)\n",
               dp_wirelength,
               ((dp_wirelength - orig_wirelength) / orig_wirelength) * 100.0);
        total_detail_time = time(NULL) - detail_time_start;
        total_time += total_detail_time;
    }

    cout << endl << endl << "////////////////////" << endl;
    if (placement.plname() != "") {
        cout << "Benchmark: " << placement.plname() << endl;
    } else {
        cout << "Benchmark: " << placement.name() << endl;
    }
    if (param.bRunGlobal) {
        printf("\nGlobal HPWL: %.0f   Time: %6.1f sec (%.1f min)\n",
               gp_wirelength,
               static_cast<double>(total_global_time),
               static_cast<double>(total_global_time) / 60.0);
    }
    if (param.bRunLegal) {
        printf(" Legal HPWL: %.0f   Time: %6.1f sec (%.1f min)\n",
               lg_wirelength,
               static_cast<double>(total_legal_time),
               static_cast<double>(total_legal_time) / 60.0);
    }
    if (param.bRunDetail && bLegal) {
        printf("Detail HPWL: %.0f   Time: %6.1f sec (%.1f min)\n",
               dp_wirelength,
               static_cast<double>(total_detail_time),
               static_cast<double>(total_detail_time) / 60.0);
    }
    printf(" ===================================================================\n");
    printf("       HPWL: %.0f   Time: %6.1f sec (%.1f min)\n",
           placement.computeHpwl(),
           static_cast<double>(total_time),
           static_cast<double>(total_time) / 60.0);

    return 0;
}
