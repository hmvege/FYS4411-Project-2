#include "wavefunctions.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

WaveFunctions::WaveFunctions(int new_nParticles, int new_nDimensions)
{
    nParticles = new_nParticles;
    nDimensions = new_nDimensions;
}

double WaveFunctions::calculate(double ** positions)
{
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::calculate" << endl;
    return 1.0; // Add stuff here?
}

double WaveFunctions::localEnergy(double ** positions)
{
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::localEnergy" << endl;
    return 1.0; // Add more stuff here as well?
}

void WaveFunctions::quantumForce(double **positions, double **F, int k)
{
    cout << "If you're seeing this, you are doing it wrong in WaveFunctions::quantumForce" << endl;
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
