#include "twoelectronjastrov.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <mpi.h>

using std::cout;
using std::endl;

twoElectronJastrov::twoElectronJastrov(int new_nParticles,
                                       int new_nDimensions,
                                       int new_numprocs,
                                       int new_processRank,
                                       double new_omega,
                                       double new_alpha,
                                       double new_a,
                                       double new_beta) : WaveFunctions(new_nParticles, new_nDimensions, new_numprocs, new_processRank)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    setOmega(new_omega);
    set_a(new_a);
    setAlpha(new_alpha);
    setBeta(new_beta);
}

void twoElectronJastrov::initializeWFSampling(double **r)
{
    /*
     *  Not used by the 2 electron WF with Jastrov sampling. Needed to be able to generalize the n-electron case.
     */
}

double twoElectronJastrov::initializeWaveFunction(double **r)
{
    /*
     * Calculates the first wavefunction without anything extra. Needed to be able to generalize the n-electron case.
     * Arguments:
     *  r   : particle positions
     */
    return calculate(r,0);
}


double twoElectronJastrov::calculate(double ** r, int k)
{
    /*
     * Calculates the wavefunction with a Jastrov factor.
     * Arguments:
     *  r   : particle positions
     *  k   : index of particle being moved
     */
    double r1 = r[0][0]*r[0][0] + r[0][1]*r[0][1]; // x1^2 + y1^2
    double r2 = r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x2^2 + y2^2
    return exp( - 0.5*omega*alpha*(r1+r2) + a/(1.0/r_ij(r[0],r[1]) + beta) );
}

void twoElectronJastrov::localEnergy(double ** r, double &ETotal, double &EKinetic, double &EPotential)
{
    /*
     * Function for calculating the two electron local energy with Jastrov factor and Coulomb interaction.
     * Arguments:
     *  r   : particle positions
     */
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    double r12 = r_ij(r[0],r[1]);
    double r12Beta = 1 + beta*r12;
    double r12BetaSquared = r12Beta*r12Beta;
    double kineticEnergy = - 0.5*( alpha*alpha*omega*omega*(rr) - 4*alpha*omega + 2*a/r12BetaSquared*( a/r12BetaSquared - omega*alpha*r12 + 1/r12 - 2*beta/r12Beta ));
    double potentialEnergy = 0.5 * omega*omega*rr;
    if (coulombInteraction)
    {
        potentialEnergy += coulomb(r);
    }
    EKinetic = kineticEnergy;
    EPotential = potentialEnergy;
    ETotal = kineticEnergy + potentialEnergy;
}

void twoElectronJastrov::quantumForce(double **r, double **F, int k)
{
    /*
     * Function for calculating the two electron quantum force.
     * Arguments:
     *  r   : particle positions
     *  F   : particle forces
     *  k   : particle being moved
     */
    double r12 = r_ij(r[0],r[1]);
    double r12Overr12r12Beta = a/(r12*(1 + r12*beta)*(1 + r12*beta));
    if (k==0)
    {
        F[k][0] = (- omega*alpha*r[k][0] + (r[k][0] - r[k][1])*r12Overr12r12Beta)*2.0; // Hardcoded to 2 electron case
        F[k][1] = (- omega*alpha*r[k][1] + (r[k][0] - r[k][1])*r12Overr12r12Beta)*2.0;
    }
    else
    {
        F[k][0] = (- omega*alpha*r[k][0] - (r[k][0] - r[k][1])*r12Overr12r12Beta)*2.0; // Hardcoded to 2 electron case
        F[k][1] = (- omega*alpha*r[k][1] - (r[k][0] - r[k][1])*r12Overr12r12Beta)*2.0;
    }
}

void twoElectronJastrov::steepestDescent(double &ESum, int NCycles)
{
    /*
     * Using steepest descent to update the variational parameters of the wavefunction.
     * Arguments:
     *  ESum        : Local energy sum from running N MC cycles
     *  NCyclces    : Number of MC cycles used
     */
    SDStatistics(NCycles);
    double alphaDerivative = 2*(dPsiEAlphaSum - dPsiAlphaSum*ESum);
    double betaDerivative = 2*(dPsiEBetaSum - dPsiBetaSum*ESum);
    alpha -= SDStepLength*alphaDerivative; // Updating alpha and beta
    beta  -= SDStepLength*betaDerivative;
    dPsiEAlphaSum = 0; // Shouldn't these be reset?
    dPsiAlphaSum = 0;
    dPsiBetaSum = 0;
    dPsiEBetaSum = 0;
}

void twoElectronJastrov::sampleSD(double **r, double &E)
{
    /*
     * Sampling used by the steepest descent algorithm.
     * Arguments:
     *  r   : particle positions
     *  E   : local energy of current positions
     */
    double rr       = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    double r12      = r_ij(r[0],r[1]);
    double r12beta  = 1 + beta*r12;
    dPsiAlpha       = - 0.5*omega*rr;                   // Derivative of WF w.r.t. alpha
    dPsiBeta        = - a*r12*r12/(r12beta*r12beta);    // Derivative of WF w.r.t. beta

    dPsiAlphaSum    += dPsiAlpha;
    dPsiEAlphaSum   += dPsiAlpha*E;
    dPsiBetaSum     += dPsiBeta;
    dPsiEBetaSum    += dPsiBeta*E;
}

void twoElectronJastrov::SDStatistics(int NCycles)
{
    /*
     * Function for retrieving Steepest Descent statistics.
     * Arguments:
     *  NCycles     : Monte Carlo cycles
     */
    dPsiBetaSum     /= double(NCycles);
    dPsiAlphaSum    /= double(NCycles);
    dPsiEAlphaSum   /= double(NCycles);
    dPsiEBetaSum    /= double(NCycles);
}

void twoElectronJastrov::finalizeSD()
{
    double tempAlpha = 0;
    double tempBeta = 0;
    MPI_Reduce(&alpha, &tempAlpha, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&beta, &tempBeta, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    alpha = tempAlpha/double(numprocs);
    beta = tempBeta/double(numprocs);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void twoElectronJastrov::printVariationalParameters(int i)
{
    /*
     * Temporary function for printing the variational parameters used.
     * Arguments:
     *  i   : current steepest descent iteration
     */
    cout << "i = " << std::setw(5) << i << " Alpha = " << std::setw(10) << alpha << " Beta = " << std::setw(10) << beta << endl;
}

void twoElectronJastrov::printUpdatedVariationalParameters()
{
    /*
     * Prints updated variational parameters from steepest descent.
     */
    cout << "Updated variational parameters: " << endl;
    cout << "Alpha = " << std::setw(10) << alpha << endl;
    cout << "Beta = " << std::setw(10) << beta << endl;
}

std::string twoElectronJastrov::getParameterString()
{
    /*
     * Returns string to be used in filename.
     */
    return "_omega" + std::to_string(omega) + "_alpha" + std::to_string(alpha) + "_beta" + std::to_string(beta);
}

void twoElectronJastrov::updateWF()
{
    /*
     * Not needed by the two electron jastrow hardcoded case.
     */
}

void twoElectronJastrov::revert()
{
    /*
     * Not needed by the two electron jastrow hardcoded case.
     */
}

void twoElectronJastrov::reset()
{
    /*
     * Not needed by the two electron jastrow hardcoded case.
     */
}
