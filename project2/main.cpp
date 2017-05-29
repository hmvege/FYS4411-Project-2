#include <iostream>
#include <cmath>
#include <ctime>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"
#include "wavefunctions/twoelectronjastrov.h"
#include "wavefunctions/n_electron_wf/nelectron.h"
#include "samplers/metropolissampler.h"
#include "samplers/uniformsampling.h"
#include "samplers/importancesampler.h"
#include "functions.h"

using namespace std;

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength, double seed, bool impSampling, bool coulomb);
void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions, double omega, double alpha, double a, double beta, double D, double deltat, double seed, double SDStepLength, bool runImpSampling, bool coulomb);
void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed,
                   double SDStepLength, bool impSampling, bool coulomb, bool jastrow);

/*
 * TODO GENERAL:
 * [ ] Switch to new r_if(r1,r2) function in all classes
 * [ ] Clean up unused variables(e.g. n_VarParams) <- CHECK ALL WARNINGS
 * [ ] Add destructors to where it is appropriate
 * [ ] Make the variable a for the two-electron cases non-global
 * [ ] Clean up main to a simpler interface
 * [ ] Add data sampling to file for every MCCount % 10000
 * [ ] Parallelize code
 */


int main()
{
    // Constants
    unsigned int MCCycles   = 1e6;
    unsigned int optCycles  = 1e5;
    int maxSDIterations     = 0;
    int nParticles          = 2;
    int nDimensions         = 2;
    double omega            = 1.0;
    double alpha            = 1.0;//0.988559;
//    double alpha            = 0.762647; // No jastrow, 2 electrons
//    double alpha            = 0.6; // No jastrow, 6 electrons
    double a                = 1.0;
    double beta             = 0.4;//0.398665;
    double D                = 0.5; // equals 0.5 in atomic units
    double deltat           = 0.001; // should be either 0.01-0.001
    double SDStepLength     = 0.01; // Steepest descent step length
    double seed             = -1;//std::time(nullptr);
    bool importanceSampling = true;
    bool coulombInteraction = true;
    bool jastrowFactor      = true;

    clock_t programStart, programEnd;
    programStart = clock();

//    run2Electron(MCCycles, nParticles, nDimensions, omega, alpha, 1.31, seed, importanceSampling, coulombInteraction);
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, a, beta, D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction);
    runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, beta, D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction, jastrowFactor);

    programEnd = clock();
    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

void printRunInfo(int MCCycles, int maxSDIterations, int nParticles, int importanceSampling, int coulombInteraction, bool jastrov)
{
    /*
     * Small function for printing out the configuration a run
     */
    for (int i = 0; i < 1e2; i++) { cout << "="; } cout << endl; // Printing a line
    cout << "N electrons                " << nParticles << endl;
    cout << "MC-cycles:                 " << MCCycles << endl;
    cout << "Jastrov                    ";
    if (jastrov) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Steepest descent:          ";
    if (maxSDIterations > 0) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Importance sampling:       ";
    if (importanceSampling) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "Coulomb interaction:       ";
    if (coulombInteraction) { cout << "TRUE" << endl; } else { cout << "FALSE" << endl; }
    cout << "RESULTS:" << endl;
}

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength, double seed, bool impSampling, bool coulomb)
{
    /*
     * Function for running the two electron case. CAN BE CONFIGURED INTO A UNIT TEST LATER
     */
    printRunInfo(MCCycles, 0, nParticles, impSampling, coulomb, false);
    twoElectronPlain WF_2Electron(nParticles, nDimensions, 1, omega, alpha);
    WF_2Electron.setCoulombInteraction(coulomb);
    VMC VMC_2Electron(nParticles, nDimensions);

    if (impSampling)
    {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Electron);
        importanceSampling.initializeSampling(0.001, seed, 0.5);
//        importanceSampling.setWaveFunction(&WF_2Electron);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Electron);
        uniformSampling.initialize(stepLength, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.runVMC(MCCycles, 0, 0);
    VMC_2Electron.printResults();
}

void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool impSampling, bool coulomb)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, true);
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, 2, omega, alpha, a, beta);
    WF_2Jastrov.setSDStepLength(SDStepLength);
    VMC VMC_2Electron(nParticles, nDimensions);

    WF_2Jastrov.setCoulombInteraction(coulomb);

//    // Temporary fix for the sampling of steepest descent metropolis part ========================
//    UniformSampling SDR(nParticles, nDimensions, &WF_2Jastrov);
//    SDR.initialize(1.14, seed-1);
//    VMC_2Electron.SDR = &SDR;
//    // ===========================================================================================

    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Jastrov);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Jastrov);
        uniformSampling.initialize(1.14, seed); // CHECK STEP!
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }

    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
    VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD);
    VMC_2Electron.printResults();
}

void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed,
                   double SDStepLength, bool impSampling, bool coulomb, bool jastrow)
{
    printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, jastrow);
    NElectron WF_NElectron(nParticles, nDimensions, 2, omega, alpha, beta);
    WF_NElectron.setSDStepLength(SDStepLength);
    VMC VMC_NElectron(nParticles, nDimensions);
    WF_NElectron.setCoulombInteraction(coulomb);
    WF_NElectron.setJastrow(jastrow);
    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_NElectron);
        importanceSampling.initializeSampling(deltat, seed, D);        
        VMC_NElectron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_NElectron);
        uniformSampling.initialize(1.14, seed); // CHECK STEP!
        VMC_NElectron.setMetropolisSampler(&uniformSampling);
    }
    VMC_NElectron.setWaveFunction(&WF_NElectron);
    VMC_NElectron.runVMC(MCCycles,optCycles,maxNSD);
    VMC_NElectron.printResults();
}
