#include <iostream>
#include <cmath>
#include "vmc.h"
#include "wavefunctions/twoelectronplain.h"

using namespace std;

double trialWF2Electron(double **rPos);
double localEnergyTwoElectron(double **rPos);

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
    twoElectronPlain WF_2Electron(nParticles, nDimensions, omega, alpha, C);

    VMC VMC_2Electron;
    VMC_2Electron.setNParticles(nParticles);
    VMC_2Electron.setNDimensions(nDimensions);
    VMC_2Electron.setWaveFunction(&WF_2Electron);
    VMC_2Electron.runVMC(MCCycles);

    programEnd = clock();
    cout << "Program complete. Time used: " << ((programEnd - programStart)/((double)CLOCKS_PER_SEC)) << endl;
    return 0;
}

double trialWF2Electron(double ** rPos)
{
    /*
     * This should be made a into a generalized WF class later
     */
    double omega = 1.0;
    double a = 1.0;
    double alpha = 1.0;
    double beta = 1.0;
    double C = 1.0;
    double x1 = rPos[0][0];
    double y1 = rPos[0][1];
    double x2 = rPos[1][0];
    double y2 = rPos[1][1];

    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);

    // With Jastrov factor
//    return C*exp( - 0.5*omega*alpha*(r12Squared) + a*r12/(1.0 + beta*r12) );

    // Without Jastrov factor
    return C*exp( - 0.5*omega*alpha*(r12Squared));
}

double localEnergyTwoElectron(double ** rPos)
{
    double omega = 1.0;
    double a = 1.0;
    double alpha = 1.0;
    double beta = 1.0;
    double x1 = rPos[0][0];
    double y1 = rPos[0][1];
    double x2 = rPos[1][0];
    double y2 = rPos[1][1];

    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);
    double r12Beta = 1 + beta*r12;
    double r12BetaSquared = r12Beta*r12Beta;

    // With Jastrov factor
//    return - 0.5*( (alpha*alpha - 1)*omega*omega*(r1*r1 + r2*r2) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta ) );

    // Without Jastrov factor. Should equal 2
    return -0.5*omega*(omega*(r12Squared) - 4) + 0.5 * omega*omega*r12Squared;
}
