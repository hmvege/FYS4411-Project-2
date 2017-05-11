#include <iostream>
#include <cmath>
#include <ctime>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"
#include "wavefunctions/twoelectronjastrov.h"
#include "samplers/metropolissampler.h"
#include "samplers/uniformsampling.h"
#include "samplers/importancesampler.h"

using namespace std;

void run2Electron(int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double C, double stepLength, double seed);
void run2eImpSampling(int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double C, double a, double beta, double D, double deltat, double seed);

int main()
{
    // Constants
    int MCCycles    = 1e6;
    int nParticles  = 2;
    int nDimensions = 2;

    clock_t programStart, programEnd;
    programStart = clock();

    // TASK C-F CONSTANTS
    double omega    = 1.0;
    double alpha    = 1.0; // 0.97
    double C        = 1.0;
    double a        = 1.0;
    double beta     = 1.0; // 0.4
    double D = 0.5; // equals 0.5 in atomic units
    double deltat = 0.001; // should be either 0.01-0.001
    double seed = std::time(nullptr);

    run2Electron(MCCycles, nParticles, nDimensions, omega, alpha, C, 1.05, seed);
    run2eImpSampling(MCCycles, nParticles, nDimensions, omega, alpha, C, a, beta, D, deltat, seed);

    programEnd = clock();
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

void run2Electron(int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double C, double stepLength, double seed)
{
    /*
     * Function for running the two electron case. CAN BE CONFIGURED INTO A UNIT TEST LATER
     */
    cout << "============ Running for 2 electron case with no Jastrov factor and uniform sampling ============" << endl;
    twoElectronPlain WF_2Electron(nParticles, nDimensions, omega, alpha, C);
    UniformSampling uniformSampling(nParticles, nDimensions);
    uniformSampling.initialize(stepLength, seed);

    VMC VMC_2Electron(nParticles, nDimensions);
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.setMetropolisSampler(&uniformSampling);
    VMC_2Electron.runVMC(MCCycles);
    VMC_2Electron.getStatistics();
}

void run2eImpSampling(int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double C, double a, double beta, double D, double deltat, double seed)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    cout << "============ Running for 2 electron case with Jastrov factor and importancesampling =============" << endl;
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, omega, alpha, C, a, beta);

    ImportanceSampler importanceSampling(nParticles, nDimensions);
    importanceSampling.initializeSampling(deltat, seed, D);
    importanceSampling.setWaveFunction(&WF_2Jastrov);

//    UniformSampling uniformSampling;
//    uniformSampling.initialize(1.14, seed); // CHECK STEP!

    VMC VMC_2Electron(nParticles, nDimensions);
    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
    VMC_2Electron.setMetropolisSampler(&importanceSampling);
//    VMC_2Electron.setMetropolisSampler(&uniformSampling);
    VMC_2Electron.runVMC(MCCycles);
    VMC_2Electron.getStatistics();

}
