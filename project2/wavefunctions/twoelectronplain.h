#ifndef TWOELECTRONPLAIN_H
#define TWOELECTRONPLAIN_H

#include "wavefunctions.h"
#include <iostream>
#include <iomanip>

class twoElectronPlain : public WaveFunctions
{
private:
    double omega;
    double alpha;
public:
    using WaveFunctions::localEnergy;

    twoElectronPlain(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank, double new_omega, double new_alpha);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    virtual double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
    std::string getParameterString();
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    // Printers
    void printVariationalParameters(int i);
};

#endif // TWOELECTRONPLAIN_H
