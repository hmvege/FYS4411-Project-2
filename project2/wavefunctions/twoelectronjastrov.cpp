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
    /*
     * Function for calculating the two electron local energy
     */
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
    /*
     * Function for calculating the two electron quantum force
     */
    double r12 = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
    double r12r12Beta = r12*(1 + beta*r12);
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

void twoElectronJastrov::steepestDescent(double ESum, int NCycles)
{
    /*
     * Should update the variational parameters of the wavefunctio.
     */
    double epsilon = 0.001; // CHANGE LATER!!
    // Statistics
    SDStatistics(NCycles);
    ESum /= double(nParticles*NCycles);
    // Getting derivatives
    double alphaDerivative = 2*(dPsiEAlphaSum - dPsiAlphaSum*ESum);
    double betaDerivative = 2*(dPsiEBetaSum - dPsiBetaSum*ESum);
    // Updating alpha and beta
    alpha -= epsilon*alphaDerivative;
    beta -= epsilon*betaDerivative;
}

void twoElectronJastrov::sampleSD(double **r, double E)
{
    /*
     * Sampling used by the steepest descent algorithm
     */
    double rr       = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // r_1^2 + r_2^2
    double r12      = sqrt((r[0][0]-r[1][0])*(r[0][0]-r[1][0]) + (r[0][1]-r[1][1])*(r[0][1]-r[1][1])); // sqrt((x1-x2)^2 + (y1-y2)^2)
    double r12beta  = 1 + beta*r12;
    dPsiAlpha       = - 0.5*omega*rr;
    dPsiAlphaSum    += dPsiAlpha;
    dPsiEAlphaSum   += dPsiAlpha*E;
    dPsiBeta        = - a*r12*r12/(r12beta*r12beta);
    dPsiBetaSum     += dPsiBeta;
    dPsiEBetaSum    += dPsiBeta*E;
}

void twoElectronJastrov::SDStatistics(int NCycles)
{
    /*
     * Function for retrieving Steepest Descent statistics. Arguments:
     * r        : positions
     * NCycles  : Monte Carlo cycles
     */
    dPsiBetaSum     /= double(nParticles*NCycles);
    dPsiAlphaSum    /= double(nParticles*NCycles);
    dPsiEAlphaSum   /= double(nParticles*NCycles);
    dPsiEBetaSum    /= double(nParticles*NCycles);
}

void twoElectronJastrov::printVariationalParameters()
{
    /*
     * Temporary function for printing the variational parameters used.
     */
    cout << "Alpha = " << alpha << " Beta = " << beta << endl;
}
