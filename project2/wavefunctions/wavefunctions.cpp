#include "wavefunctions.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

WaveFunctions::WaveFunctions(int new_nParticles, int new_nDimensions, int new_nVarParams)
{
    cout << "TODO: make a general relative distance function" << endl;
    nParticles = new_nParticles;
    nDimensions = new_nDimensions;
    nVarParams = new_nVarParams;
}

double WaveFunctions::calculate(double ** r)
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
                r12abs += (r[i][k] - r[j][k])*(r[i][k] - r[j][k]);
            }
            coulombInteraction += 1/sqrt(r12abs);
        }
    }
    return coulombInteraction;
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
