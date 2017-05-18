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

double twoElectronPlain::calculate(double ** r)
{
    /*
     * Calculates the wavefunction
     */
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    return exp( - 0.5*omega*alpha*(rr)); // No Jastrov-factor
}

double twoElectronPlain::localEnergy(double ** r)
{
    /*
     * Calculates local energy of the two electron case without any Jastrov factor and no Coulomb interaction
     */
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    return -0.5*omega*(omega*(rr) - 4) + 0.5 * omega*omega*(rr); // No Jastrov-factor, no Coulomb interaction
}

void twoElectronPlain::printVariationalParameters()
{
    cout << "Alpha = " << alpha << endl;
}

