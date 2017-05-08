#include "twoelectronplain.h"
#include <cmath>

twoElectronPlain::twoElectronPlain(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha, double new_C)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    nParticles  = new_nParticles;
    nDimensions = new_nDimensions;
    omega       = new_omega;
    alpha       = new_alpha;
    C           = new_C;
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
    double r12Squared = r1*r1 + r2*r2;
    return C*exp( - 0.5*omega*alpha*(r12Squared)); // No Jastrov-factor
}

double twoElectronPlain::localEnergy(double ** positions)
{
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    return -0.5*omega*(omega*(r12Squared) - 4) + 0.5 * omega*omega*r12Squared; // No Jastrov-factor, no Coulomb interaction
}

double **twoElectronPlain::quantumForce(double **positions)
{
    return positions; // Really does nothing by default
}