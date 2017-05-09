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

double twoElectronJastrov::localEnergy(double ** positions)
{
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);
    double r12Beta = 1 + beta*r12;
    double r12BetaSquared = r12Beta*r12Beta;

//    double coulombInteraction = 1/sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
    // General method for getting the coulomb interaction value
    double coulombInteraction = 0.0;
    double r12abs = 0.0;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < i; j++)
        {
            r12abs = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r12abs += sqrt( (positions[i][k] - positions[j][k])*(positions[i][k] - positions[j][k]) );
            }
            coulombInteraction += 1/r12abs;
        }
    }
//    cout << "Calculating local energy: remember to check that the coloumb interaction is correct" << endl;
    // With Jastrov factor and Coulomb interaction
    return - 0.5*( (alpha*alpha - 1)*omega*omega*(r1*r1 + r2*r2) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta )) + coulombInteraction;
}

double *twoElectronJastrov::quantumForce(double **positions, int k)
{
    double x1 = positions[0][0];
    double y1 = positions[0][1];
    double x2 = positions[1][0];
    double y2 = positions[1][1];
//    cout << "inside jastrov quantum force"<< endl;
    double r1 = x1*x1 + y1*y1;
    double r2 = x2*x2 + y2*y2;
    double r12Squared = r1*r1 + r2*r2;
    double r12 = sqrt(r12Squared);
    double r12Beta = 1 + beta*r12;
    double r12r12Beta = r12*(r12Beta);
    double * F = new double [nDimensions]; // Force for k particle
    for (int i = 0; i < nDimensions; i++)
    {
        F[i] = 0;
    }
    for (int i = 0; i < nDimensions; i++)
    {
        F[i] = - 2*omega*alpha*positions[k][i];
        if ((i + 1) % 2 == 0) // If index == 0
        {
            F[i] -= a*(positions[k][0] - positions[k][1])/r12r12Beta;
        }
        else // If index == 1
        {
            F[i] += a*(positions[k][0] - positions[k][1])/r12r12Beta;
        }
    }
    return F; // FIX THIS!!
}

