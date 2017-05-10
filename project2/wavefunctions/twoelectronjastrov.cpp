#include "twoelectronjastrov.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

twoElectronJastrov::twoElectronJastrov(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha, double new_C, double new_a, double new_beta)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    nParticles  = new_nParticles;
    nDimensions = new_nDimensions;
    omega       = new_omega;
    a           = new_a;
    alpha       = new_alpha;
    beta        = new_beta;
    C           = new_C;
}

double twoElectronJastrov::calculate(double ** r)
{
    /*
     * Calculates the wavefunction with a Jastrov factor.
     */
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);
    return C*exp( - 0.5*omega*alpha*(r12Squared) + a*r12/(1.0 + beta*r12) );
}

double twoElectronJastrov::coulomb(double ** r)
{
    // General method for getting the coulomb interaction value
    double coulombInteraction = 0.0; // PUT THIS INTO CLASS
    double r12abs;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < i; j++)
        {
            r12abs = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r12abs += sqrt( (r[i][k] - r[j][k])*(r[i][k] - r[j][k]) );
            }
            coulombInteraction += 1/r12abs;
        }
    }
    return coulombInteraction;
}

double twoElectronJastrov::localEnergy(double ** r)
{
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    double r12 = sqrt(r1*r1 + r2*r2);
    double r12Beta = 1 + beta*r12;
    double r12BetaSquared = r12Beta*r12Beta;

    double coulombInteraction = 1/(sqrt((r[1][0]-r[0][0])*(r[1][0]-r[0][0])) + sqrt((r[1][1]-r[0][1])*(r[1][1]-r[0][1]))); // Hardcoded 2 electron case

    // With Jastrov factor and Coulomb interaction
    return - 0.5*( (alpha*alpha - 1)*omega*omega*(r1*r1 + r2*r2) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta )) + coulombInteraction;
}

double *twoElectronJastrov::quantumForce(double **r, int k)
{
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    double r12 = sqrt(r1*r1 + r2*r2);
    double r12r12Beta = r12*(1 + beta*r12);
    double * F = new double [nDimensions]; // Force for k particle
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
    F[0] = (- omega*alpha*r[k][0] + a*(r[k][0] - r[k][1])/r12r12Beta)*2.0; // Hardcoded to 2 electron case
    F[1] = (- omega*alpha*r[k][1] - a*(r[k][0] - r[k][1])/r12r12Beta)*2.0;
    return F; // FIX THIS!!
}

