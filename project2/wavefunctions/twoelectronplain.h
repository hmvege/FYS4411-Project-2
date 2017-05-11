#ifndef TWOELECTRONPLAIN_H
#define TWOELECTRONPLAIN_H

#include "wavefunctions.h"

class twoElectronPlain : public WaveFunctions
{
private:
    double omega;
    double alpha;
    double C;
public:
//    twoElectronPlain() { }
    twoElectronPlain(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha, double new_C);

    double calculate(double **positions);
    double localEnergy(double **positions);
//    double *quantumForce(double **positions);
};

#endif // TWOELECTRONPLAIN_H
