#ifndef TWOELECTRONPLAIN_H
#define TWOELECTRONPLAIN_H

#include "wavefunctions.h"

class twoElectronPlain : public WaveFunctions
{
private:
    double omega;
    double alpha;
public:
    twoElectronPlain(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha);
    double calculate(double **positions);
    double localEnergy(double **positions);
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
};

#endif // TWOELECTRONPLAIN_H
