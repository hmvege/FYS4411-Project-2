#include "twoelectronplain.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

twoElectronPlain::twoElectronPlain(int new_nParticles,
                                   int new_nDimensions,
                                   int new_nVarParams,
                                   double new_omega,
                                   double new_alpha) : WaveFunctions(new_nParticles, new_nDimensions, new_nVarParams)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    nParticles  = new_nParticles;
    nDimensions = new_nDimensions;
    omega       = new_omega;
    alpha       = new_alpha;
}

double twoElectronPlain::calculate(double ** positions)
{
    /*
     * Calculates the wavefunction
     */
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
//    double r12Squared = r1 + r2*r2;
    return exp( - 0.5*omega*alpha*(r1 + r2)); // No Jastrov-factor
}

double twoElectronPlain::localEnergy(double ** positions)
{
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
//    double r12Squared = r1*r1 + r2*r2;
    return -0.5*omega*(omega*(r1 + r2) - 4) + 0.5 * omega*omega*(r1+r2); // No Jastrov-factor, no Coulomb interaction
}

void twoElectronPlain::printVariationalParameters()
{
    cout << "Alpha = " << alpha << endl;
}

void twoElectronPlain::steepestDescent(double **r)
{
    double epsilon = 0.001; // CHANGE LATER!!
    double alphaDerivative = omega*(2-omega*alpha*(r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]));
    alpha -= epsilon*alphaDerivative;
}

//double *twoElectronPlain::quantumForce(double **positions)
//{
//    cout << "twoelectronplain.cpp: quantumForce not implemented" << endl;
//    return nullptr; // Really does nothing by default
//}
