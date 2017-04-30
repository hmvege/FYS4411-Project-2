#include "twoelectronjastrov.h"
#include <cmath>

twoElectronJastrov::twoElectronJastrov(double omega_, double a_, double alpha_, double beta_, double C_)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    omega   = omega_;
    a       = a_;
    alpha   = alpha_;
    beta    = beta_;
    C       = C_;
}

double twoElectronJastrov::calculate(double ** positions)
{
    /*
     * Calculates the wavefunction with a Jastrov factor.
     */
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);
    return C*exp( - 0.5*omega*alpha*(r12Squared) + a*r12/(1.0 + beta*r12) );
}
