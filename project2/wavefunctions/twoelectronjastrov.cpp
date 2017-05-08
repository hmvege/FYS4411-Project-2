#include "twoelectronjastrov.h"
#include <cmath>

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

double twoElectronJastrov::localEnergy(double ** rPos)
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
                r12abs += sqrt( (rPos[i][k] - rPos[j][k])*(rPos[i][k] - rPos[j][k]) );
            }
            coulombInteraction += 1/r12abs;
        }
    }

    // With Jastrov factor and Coulomb interaction
    return - 0.5*( (alpha*alpha - 1)*omega*omega*(r1*r1 + r2*r2) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta )) + coulombInteraction;
}

double **twoElectronJastrov::quantumForce(double **positions)
{
    double ** F = new double * [nParticles];
    for (int i = 0; i < nParticles; i++) // General method that can be used for later
    {
        F[i] = new double [nDimensions];
        for (int j = 0; j < nDimensions; j++)
        {
            F[i][j] = 0;
        }
    }
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < nDimensions; j++)
        {
            for (int k = 0; k < nParticles; k++)
            {
                F[i][j] += positions[i][k];
            }
            F[i][j] *= -2*omega*alpha;
        }
    }
    return F; // FIX THIS!!
}

