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
//    void initialize(double **r, double &WF);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
//    void steepestDescent(double E, int NCycles);
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    // Printers
    void printVariationalParameters();
};

#endif // TWOELECTRONPLAIN_H
