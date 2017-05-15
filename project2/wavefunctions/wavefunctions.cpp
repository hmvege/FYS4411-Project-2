#include "wavefunctions.h"
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
