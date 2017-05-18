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

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength, double seed);
void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions, double omega, double alpha, double a, double beta, double D, double deltat, double seed, double SDStepLength, bool runImpSampling);

int main()
{
    // Constants
    unsigned int MCCycles   = 1e6;
    unsigned int optCycles  = 1e5;
    int maxSDIterations     = 0;
    int nParticles          = 2;
    int nDimensions         = 2;

    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    clock_t programStart, programEnd;
    programStart = clock();

    // TASK C-F CONSTANTS
    double omega            = 1.0;
    double alpha            = 0.988559;
    double a                = 1.0;
    double beta             = 0.398665;
    double D                = 0.5; // equals 0.5 in atomic units
    double deltat           = 0.001; // should be either 0.01-0.001
    double SDStepLength     = 0.01; // Steepest descent step length
    double seed             = std::time(nullptr);

//    run2Electron(MCCycles, nParticles, nDimensions, omega, alpha, 1.31, seed);
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, a, beta, D, deltat, seed, SDStepLength, false);
    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, a, beta, D, deltat, seed, SDStepLength, true);

    programEnd = clock();
    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength, double seed)
{
    /*
     * Function for running the two electron case. CAN BE CONFIGURED INTO A UNIT TEST LATER
     */
    cout << "============== Running for 2 electron case with no Jastrov factor and uniform sampling =============" << endl;
    twoElectronPlain WF_2Electron(nParticles, nDimensions, 1, omega, alpha);
    UniformSampling uniformSampling(nParticles, nDimensions);
    uniformSampling.initialize(stepLength, seed);

    VMC VMC_2Electron(nParticles, nDimensions);
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.setMetropolisSampler(&uniformSampling);
    VMC_2Electron.runVMC(MCCycles, 0, 0);
    VMC_2Electron.printResults();
}

void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool runImpSampling)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, 2, omega, alpha, a, beta);
    WF_2Jastrov.setSDStepLength(SDStepLength);
    VMC VMC_2Electron(nParticles, nDimensions);

    // Temporary fix for the sampling of steepest descent metropolis part ========================
    UniformSampling SDR(nParticles, nDimensions);
    SDR.initialize(1.14, seed-1);
    VMC_2Electron.SDR = &SDR;
    // ===========================================================================================

    if (runImpSampling) {
        cout << "============== Running for 2 electron case with Jastrov factor and importancesampling ==============" << endl;
        ImportanceSampler importanceSampling(nParticles, nDimensions);
        importanceSampling.initializeSampling(deltat, seed, D);
        importanceSampling.setWaveFunction(&WF_2Jastrov);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        cout << "=============== Running for 2 electron case with Jastrov factor and uniform sampling ===============" << endl;
        UniformSampling uniformSampling(nParticles, nDimensions);
        uniformSampling.initialize(1.14, seed); // CHECK STEP!
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }

    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
    VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD);
    VMC_2Electron.printResults();

}
