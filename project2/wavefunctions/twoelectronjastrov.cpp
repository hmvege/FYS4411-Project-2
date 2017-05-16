#include "twoelectronjastrov.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

twoElectronJastrov::twoElectronJastrov(int new_nParticles,
                                       int new_nDimensions,
                                       int new_nVarParams,
                                       double new_omega,
                                       double new_alpha,
                                       double new_a,
                                       double new_beta) : WaveFunctions(new_nParticles, new_nDimensions, new_nVarParams)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
//    nParticles  = new_nParticles;
//    nDimensions = new_nDimensions;
    omega       = new_omega;
    a           = new_a;
    alpha       = new_alpha;
    beta        = new_beta;
}

double twoElectronJastrov::calculate(double ** r)
{
    /*
     * Calculates the wavefunction with a Jastrov factor.
     */
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    double r12 = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((r1x-r2x)^2 + (r1y-r2y)^2)
    return exp( - 0.5*omega*alpha*(r1+r2) + a*r12/(1.0 + beta*r12) );
}

double twoElectronJastrov::localEnergy(double ** r)
{
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    double r12 = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
    double r12Beta = 1 + beta*r12;
    double r12BetaSquared = r12Beta*r12Beta;
    // With Jastrov factor and Coulomb interaction
    return - 0.5*( (alpha*alpha - 1)*omega*omega*(r1 + r2) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta )) + coulomb(r);
}

void twoElectronJastrov::quantumForce(double **r, double **F, int k)
{
    double r12 = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
    double r12r12Beta = r12*(1 + beta*r12);
//    for (int i = 0; i < nDimensions; i++)
//    {
//        F[i] = 0;
//    }
//    for (int i = 0; i < nDimensions; i++)
//    {
//        F[i] = - 2*omega*alpha*r[k][i];
//        if ((i + 1) % 2 == 0) // If index == 1, that is we get an extra minus sign
//        {
//            F[i] -= a*(r[k][0] - r[k][1])/r12r12Beta;
//        }
//        else // If index == 1
//        {
//            F[i] += a*(r[k][0] - r[k][1])/r12r12Beta;
//        }
//    }
    if (k==0) {
        F[k][0] = (- omega*alpha*r[k][0] + a*(r[k][0] - r[k][1])/r12r12Beta)*2.0; // Hardcoded to 2 electron case
        F[k][1] = (- omega*alpha*r[k][1] + a*(r[k][0] - r[k][1])/r12r12Beta)*2.0;
    }
    else
    {
        F[k][0] = (- omega*alpha*r[k][0] - a*(r[k][0] - r[k][1])/r12r12Beta)*2.0; // Hardcoded to 2 electron case
        F[k][1] = (- omega*alpha*r[k][1] - a*(r[k][0] - r[k][1])/r12r12Beta)*2.0;
    }
}

void twoElectronJastrov::steepestDescent(double **r)
{
    /*
     * Should update the variational parameters of the wavefunctio.
     */
    double epsilon = 0.001; // CHANGE LATER!!
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // r_1^2 + r_2^2
    double r12 = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
    double r12beta = (1 + beta*r12);
    double alphaDerivative = 2*omega - alpha*omega*omega*rr + omega*a*r12/(r12beta*r12beta);
    double betaDerivative = a/(r12beta*r12beta*r12beta)*( 4*r12*a/(r12beta*r12beta) - 2*omega*alpha*r12*r12 + 4 - 6*beta*r12/r12beta );
    alpha -= epsilon*alphaDerivative;
    beta -= epsilon*betaDerivative;
}

void twoElectronJastrov::printVariationalParameters()
{
    cout << "Alpha = " << alpha << " Beta = " << beta << endl;
}
