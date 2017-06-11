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
    // For steepest descent
    double dPsiAlpha = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    void SDStatistics(int NCycles);
public:
    using WaveFunctions::localEnergy;

    twoElectronPlain(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank, double new_omega, double new_alpha);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    virtual void localEnergy(double **r, double &ETotal, double &EKinetic, double &EPotential);
    void quantumForce(double **r, double **F, int k);
    void finalizeSD();
    void sampleSD(double **r, double &E);
    void steepestDescent(double &ESum, int NCycles);
    void revert();
    void updateWF();
    void reset();
    std::string getParameterString();
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    // Printers
    void printVariationalParameters(int i);
    void printUpdatedVariationalParameters();
};

#endif // TWOELECTRONPLAIN_H
