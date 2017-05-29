#include "wavefunctions.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

WaveFunctions::WaveFunctions(int new_nParticles, int new_nDimensions, int new_nVarParams)
{
    nParticles = new_nParticles;
    nDimensions = new_nDimensions;
    nVarParams = new_nVarParams;
}

void WaveFunctions::initializeWFSampling(double ** r)
{
    /*
     * By default, does nothing. Only used in the N-electron sampling case.
     *  r : position
     */
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::initializeWFClass" << endl;
    exit(1);
}

double WaveFunctions::initializeWaveFunction(double **r)
{
    /*
     * Returns the wavefunction at a given coordinate.
     *  r : position
     */
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::initializeWaveFunction" << endl;
    exit(1);
    return 1.0; // Add stuff here?

}

double WaveFunctions::calculate(double ** r, int k)
{
    /*
     * Returns the wavefunction at a given coordinate.
     * r : position
     */
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::calculate" << endl;
    exit(1);
    return 1.0; // Add stuff here?
}

double WaveFunctions::localEnergy(double ** r)
{
    /*
     * Returns the energy of the wavefunction at a given coordinate.
     * r : position
     */
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::localEnergy" << endl;
    exit(1);
    return 1.0; // Add more stuff here as well?
}

void WaveFunctions::quantumForce(double **r, double **F, int k)
{
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::quantumForce" << endl;
    exit(1);
}

double WaveFunctions::coulomb(double ** r)
{
    /*
     * General method for getting the coulomb interaction value
     */
    double coulombInteraction = 0.0; // PUT THIS INTO CLASS AS A UPPER HALF MATRIX?
    double r12abs;
    for (int i = 0; i < nParticles; i++)
    {
        for (int j = 0; j < i; j++)
        {
            r12abs = 0.0;
            for (int k = 0; k < nDimensions; k++)
            {
                r12abs += (r[i][k] - r[j][k])*(r[i][k] - r[j][k]);
            }
            coulombInteraction += 1/sqrt(r12abs);
        }
    }
    return coulombInteraction;
}

double WaveFunctions::r_ij(double *r1, double *r2)
{
    /*
     * Method for finding distance between two electrons
     */
    return sqrt((r2[0]-r1[0])*(r2[0]-r1[0]) + (r2[1]-r1[1])*(r2[1]-r1[1]));
}

void WaveFunctions::steepestDescent(double &E, int NCycles)
{
    cout << "Steepest descent not implemented for general WaveFunction module." << endl;
    exit(1);
}

void WaveFunctions::sampleSD(double **r, double &E)
{
    cout << "Steepest descent not implemented for general WaveFunction module." << endl;
    exit(1);
}

void WaveFunctions::printVariationalParameters()
{
    cout << "Bare wavefunctions class do not contain any variational parameters." << endl;
    exit(1);
}

bool WaveFunctions::SDConvergenceCriteria()
{
    cout << "Bare wavefunctions class do not Steepest Descent convergence criteria test." << endl;
    exit(1);
}

void WaveFunctions::revert(double ** r)
{
   /*
    * Function used when resetting certain class-contained variables. Needed by the n-electron case.
    */
}
