#include <iostream>
#include <cmath>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"
#include "wavefunctions/twoelectronjastrov.h"
#include "ratios/metropolisratio.h"
#include "ratios/importancesampler.h"

using namespace std;

int main()
{
    // Constants
    int MCCycles    = 1e4;
    int nParticles  = 2;
    int nDimensions = 2;

    clock_t programStart, programEnd;
    programStart = clock();

    // TASK C
    double omega    = 1.0;
    double alpha    = 1.0;
    double C        = 1.0;
    double a        = 1.0;
    double beta     = 1.0;
    twoElectronPlain WF_2Electron(nParticles, nDimensions, omega, alpha, C);
//    MetropolisRatio nonImpRatio;
    double D = 0.5; // equals 0.5 in atomic units
    double deltat = 0.001; // should be either 0.01-0.001
//    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, omega, alpha, C, a, beta);
//    ImportanceSampler ImpRatio(D, deltat, nParticles, nDimensions);
//    ImpRatio.setWaveFunction(&WF_2Jastrov);

    VMC VMC_2Electron(nParticles, nDimensions);
    VMC_2Electron.setWaveFunction(&WF_2Electron);
//    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
//    VMC_2Electron.setMetropolisRatio(&ImpRatio);
    VMC_2Electron.runVMC(MCCycles);
    VMC_2Electron.getStatistics();

    programEnd = clock();
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}
