#include "harmonicoscillator.h"

HarmonicOscillator::HarmonicOscillator(int new_nParticles, int new_nDimensions, int new_nVarParams) : WaveFunctions(new_nParticles, new_nDimensions, new_nVarParams)
{

}

//double HarmonicOscillator::calculate(double ** r)
//{
//    /*
//     * Energy of HO
//     */
//    return r[0]*r[0]*0.5 + 1.0/(8*r[0]*r[0]);
//}

//double HarmonicOscillator::localEnergy(double ** r)
//{
//    return r[0] - 1.0/(4*r[0]*r[0]*r[0]);
//}
