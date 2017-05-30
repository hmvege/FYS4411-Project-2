#include <iostream>
#include <cmath>
#include <ctime>
#include <mpi.h>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"
#include "wavefunctions/twoelectronjastrov.h"
#include "wavefunctions/n_electron_wf/nelectron.h"
#include "samplers/metropolissampler.h"
#include "samplers/uniformsampling.h"
#include "samplers/importancesampler.h"
#include "functions.h"

using namespace std;

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength,
                  double seed, bool impSampling, bool coulomb, string filename, int MCSamplingFrequency);
void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool runImpSampling, bool coulomb, string filename, int MCSamplingFrequency);
void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed, double SDStepLength,
                   bool impSampling, bool coulomb, bool jastrow, string filename, int MCSamplingFrequency);

/*
 * TODO GENERAL:
 * [x] Switch to new r_if(r1,r2) function in all classes
 * [x] Clean up unused variables(e.g. n_VarParams) <- CHECK ALL WARNINGS
 * [ ] Add data sampling to file for every MCCount % 10000
 * [ ] Parallelize code
 * [ ] Clean up main to a simpler interface for the heavy data gathering run
 * [ ] Setup blocking
 */


int main()
{
    // Constants
    unsigned int MCCycles   = 1e6;
    unsigned int optCycles  = 1e5;
    int MCSamplingFrequency = 1e5;
    int maxSDIterations     = 0;
    int nParticles          = 6;
    int nDimensions         = 2;
    double omega            = 1.0;
    double alpha            = 1.0;//0.988559; // 2 electrons
    double beta             = 0.4;//0.398665; // 2 electrons
//    double alpha            = 1.03741;    // 6 electrons
//    double beta             = 0.472513;   // 6 electrons
//    double alpha            = 1.10364;  // 12 electrons
//    double beta             = 0.468861; // 12 electrons
//    double alpha            = 0.686717; // No jastrow, 2 electrons
//    double alpha            = 0.599; // No jastrow, 6 electrons
//    double alpha            = 0.569619; // No jastrow, 12 electrons
    double a                = 1.0;
    double D                = 0.5; // equals 0.5 in atomic units
    double deltat           = 0.001; // should be either 0.01-0.001
    double SDStepLength     = 0.01; // Steepest descent step length
    double seed             = -1;//std::time(nullptr);
    bool importanceSampling = false;
    bool coulombInteraction = true;
    bool jastrowFactor      = true;

    clock_t programStart, programEnd;
    programStart = clock();

//    run2Electron(MCCycles, nParticles, nDimensions, omega, alpha, 1.31, seed, importanceSampling, coulombInteraction, "2ElectronPlain", MCSamplingFrequency);
//    run2eImpSampling(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, a, beta,D, deltat, seed, SDStepLength, importanceSampling, coulombInteraction, "2ElectronJastrov", MCSamplingFrequency);
    runNElectrons(MCCycles, optCycles, maxSDIterations, nParticles, nDimensions, omega, alpha, beta, D, deltat,seed, SDStepLength, importanceSampling, coulombInteraction, jastrowFactor, "NElectron", MCSamplingFrequency);

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

void run2Electron(unsigned int MCCycles, int nParticles, int nDimensions, double omega, double alpha, double stepLength, double seed, bool impSampling, bool coulomb, std::string filename, int MCSamplingFrequency)
{
    /*
     * Function for running the two electron case.
     */
    printRunInfo(MCCycles, 0, nParticles, impSampling, coulomb, false);
    twoElectronPlain WF_2Electron(nParticles, nDimensions, omega, alpha);
    WF_2Electron.setCoulombInteraction(coulomb);
    VMC VMC_2Electron(nParticles, nDimensions, filename);

    if (impSampling)
    {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Electron);
        importanceSampling.initializeSampling(0.001, seed, 0.5);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Electron);
        uniformSampling.initializeSampling(stepLength, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.runVMC(MCCycles, 0, 0, MCSamplingFrequency);
    VMC_2Electron.printResults();
}

void run2eImpSampling(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                      double omega, double alpha, double a, double beta, double D, double deltat, double seed,
                      double SDStepLength, bool impSampling, bool coulomb, std::string filename, int MCSamplingFrequency)
{
    /*
     * Function for running the two electron case with Jastrov factor and importance sampling.
     */
    printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, true);
    twoElectronJastrov WF_2Jastrov(nParticles, nDimensions, omega, alpha, a, beta);
    WF_2Jastrov.setSDStepLength(SDStepLength);
    VMC VMC_2Electron(nParticles, nDimensions, filename);
    WF_2Jastrov.setCoulombInteraction(coulomb);
    if (impSampling) {
        ImportanceSampler importanceSampling(nParticles, nDimensions, &WF_2Jastrov);
        importanceSampling.initializeSampling(deltat, seed, D);
        VMC_2Electron.setMetropolisSampler(&importanceSampling);
    }
    else
    {
        UniformSampling uniformSampling(nParticles, nDimensions, &WF_2Jastrov);
        uniformSampling.initializeSampling(1.14, seed);
        VMC_2Electron.setMetropolisSampler(&uniformSampling);
    }

    VMC_2Electron.setWaveFunction(&WF_2Jastrov);
    VMC_2Electron.runVMC(MCCycles, optCycles, maxNSD, MCSamplingFrequency);
    VMC_2Electron.printResults();
}

void runNElectrons(unsigned int MCCycles, unsigned int optCycles, int maxNSD, int nParticles, int nDimensions,
                   double omega, double alpha, double beta, double D, double deltat, double seed,
                   double SDStepLength, bool impSampling, bool coulomb, bool jastrow, std::string filename, int MCSamplingFrequency)
{
    printRunInfo(MCCycles, maxNSD, nParticles, impSampling, coulomb, jastrow);
    NElectron WF_NElectron(nParticles, nDimensions, omega, alpha, beta);
    WF_NElectron.setSDStepLength(SDStepLength);
    VMC VMC_NElectron(nParticles, nDimensions, filename);
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
        uniformSampling.initializeSampling(1.14, seed);
        VMC_NElectron.setMetropolisSampler(&uniformSampling);
    }
    VMC_NElectron.setWaveFunction(&WF_NElectron);
    VMC_NElectron.runVMC(MCCycles,optCycles,maxNSD, MCSamplingFrequency);
    VMC_NElectron.printResults();
}
