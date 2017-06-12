#include "twoelectronplain.h"
#include <cmath>
#include <mpi.h>
//#include <iostream>
//#include <iomanip>

using std::cout;
using std::endl;

twoElectronPlain::twoElectronPlain(int new_nParticles,
                                   int new_nDimensions,
                                   int new_numprocs,
                                   int new_processRank,
                                   double new_omega,
                                   double new_alpha) : WaveFunctions(new_nParticles, new_nDimensions, new_numprocs, new_processRank)
{
    /*
     * Class for a two-electron system. Energy should be equal to 2, and variance should be 0.
     */
    setOmega(new_omega);
    setAlpha(new_alpha);
    nParticles          = new_nParticles;
    nDimensions         = new_nDimensions;
    coulombInteraction  = false; // By default, false
}

double twoElectronPlain::initializeWaveFunction(double **r)
{
    /*
     * Calculates the first wavefunction without anything extra. Needed to be able to generalize the n-electron case.
     */
    return calculate(r,0);
}

void twoElectronPlain::initializeWFSampling(double **r)
{
    /*
     *  Not used by the plain 2 electron WF. Needed to be able to generalize the n-electron case.
     */
}

double twoElectronPlain::calculate(double ** r, int k)
{
    /*
     * Calculates the wavefunction.
     * Arguments:
     *  r   : particle positions
     *  k   : position of particle being moved
     */
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    return exp( - 0.5*omega*alpha*rr); // No Jastrov-factor
}

void twoElectronPlain::localEnergy(double ** r, double &ETotal, double &EKinetic, double &EPotential)
{
    /*
     * Calculates local energy of the two electron case without any Jastrov factor and no Coulomb interaction
     * Arguments:
     *  r   : particle positions
     */
    double rr = r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]; // x1^2 + y1^2 + x2^2 + y2^2
    double kineticEnergy = -0.5*omega*alpha*(alpha*omega*rr - 4); // with or withour alpha?
    double potentialEnergy = 0.5 * omega*omega*rr;
    if (coulombInteraction)
    {
        potentialEnergy += coulomb(r);
    }
    EKinetic = kineticEnergy;
    EPotential = potentialEnergy;
    ETotal = kineticEnergy + potentialEnergy;
}

void twoElectronPlain::quantumForce(double **r, double **F, int k)
{
    /*
     * Quantum force for the two electron case.
     * Arguments:
     *  r   : particle positions
     *  F   : quantum force vector
     *  k   : particle we are finding the force for
     */
    for (int i = 0; i < nDimensions; i++)
    {
        F[k][i] = - 2*omega*alpha*r[k][i];
    }
}

void twoElectronPlain::steepestDescent(double &ESum, int NCycles)
{
    /*
     * Using steepest descent to update the variational parameters of the wavefunction.
     * Arguments:
     *  ESum        : Local energy sum from running N MC cycles
     *  NCyclces    : Number of MC cycles used
     */
    SDStatistics(NCycles);
    double alphaDerivative = 2*(dPsiEAlphaSum - dPsiAlphaSum*ESum);
    alpha -= SDStepLength*alphaDerivative; // Updating alpha and beta
    dPsiEAlphaSum = 0; // Shouldn't these be reset?
    dPsiAlphaSum = 0;
}

void twoElectronPlain::sampleSD(double **r, double &E)
{
    /*
     * Sampling used by the steepest descent algorithm.
     * Arguments:
     *  r   : particle positions
     *  E   : local energy of current positions
     */
    dPsiAlpha       = - 0.5*omega*(r[0][0]*r[0][0] + r[0][1]*r[0][1] + r[1][0]*r[1][0] + r[1][1]*r[1][1]); // x1^2 + y1^2 + x2^2 + y2^2
    dPsiAlphaSum    += dPsiAlpha;
    dPsiEAlphaSum   += dPsiAlpha*E;
}

void twoElectronPlain::SDStatistics(int NCycles)
{
    /*
     * Function for retrieving Steepest Descent statistics.
     * Arguments:
     *  NCycles     : Monte Carlo cycles
     */
    dPsiAlphaSum    /= double(NCycles);
    dPsiEAlphaSum   /= double(NCycles);
}

void twoElectronPlain::finalizeSD()
{
    /*
     * SD for the two electron WF, more or less a carbon copy of the Jastrow case..
     */
    double tempAlpha = 0;
    MPI_Reduce(&alpha, &tempAlpha, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    alpha = tempAlpha/double(numprocs);
    MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void twoElectronPlain::printVariationalParameters(int i)
{
    /*
     * Prints variational parameters
     * Current steepest descent iteration
     */
    cout << "i = " << std::setw(5) << i << " Alpha = " << alpha << endl;
}

void twoElectronPlain::printUpdatedVariationalParameters()
{
    /*
     * Prints updated variational parameters from steepest descent.
     */
    cout << "Updated variational parameters: " << endl;
    cout << "Alpha:                     " << std::setw(10) << alpha << endl;
}

std::string twoElectronPlain::getParameterString()
{
    /*
     * Returns string to be used in filename.
     */
    return "_omega" + std::to_string(omega) + "_alpha" + std::to_string(alpha);
}

void twoElectronPlain::updateWF()
{
    /*
     * Not needed by the two electron plain hardcoded case.
     */
}

void twoElectronPlain::revert()
{
    /*
     * Not needed by the two electron plain hardcoded case.
     */
}

void twoElectronPlain::reset()
{
    /*
     * Not needed by the two electron plain hardcoded case.
     */
}
