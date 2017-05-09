#include "wavefunctions.h"

#include <iostream>

WaveFunctions::WaveFunctions(int new_nParticles, int new_nDimensions)
{
    nParticles = new_nParticles;
    nDimensions = new_nDimensions;
}

double WaveFunctions::calculate(double ** positions)

{
    return 1.0; // Add stuff here?
}

double WaveFunctions::localEnergy(double ** positions)
{
    return 1.0; // Add more stuff here as well?
}

double *WaveFunctions::quantumForce(double **positions, int k)
{
    std::cout << "If you're seeing this, you are doing it wrong" << std::endl;
    return nullptr; // Really does nothing by default
}
