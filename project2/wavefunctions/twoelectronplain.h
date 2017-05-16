#ifndef TWOELECTRONPLAIN_H
#define TWOELECTRONPLAIN_H

#include "wavefunctions.h"

class twoElectronPlain : public WaveFunctions
{
private:
    double omega;
    double alpha;
public:
    twoElectronPlain(int new_nParticles, int new_nDimensions, int new_nVarParams, double new_omega, double new_alpha);
    double calculate(double **positions);
    double localEnergy(double **positions);
    void steepestDescent(double **r);
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    // Printers
    void printVariationalParameters();
};

#endif // TWOELECTRONPLAIN_H
