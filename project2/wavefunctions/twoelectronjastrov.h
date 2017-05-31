#ifndef TWOELECTRONJASTROV_H
#define TWOELECTRONJASTROV_H

#include "wavefunctions.h"
#include <iostream>

class twoElectronJastrov : public WaveFunctions
{
private:
    double omega;
    double a;
    double alpha;
    double beta;
    // For steepest descent
    double dPsiAlpha = 0;
    double dPsiBeta = 0;
    double dPsiBetaSum = 0;
    double dPsiAlphaSum = 0;
    double dPsiEAlphaSum = 0;
    double dPsiEBetaSum = 0;
    void SDStatistics(int NCycles);
public:
    using WaveFunctions::localEnergy;
    twoElectronJastrov(int new_nParticles, int new_nDimensions, int new_numprocs, int new_processRank, double new_omega, double new_alpha, double new_a, double new_beta);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    virtual double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
    void steepestDescent(double &ESum, int NCycles);
    void sampleSD(double **r, double &E);
    std::string getParameterString();
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }
    void set_a(double new_a) { a = new_a; }
    // Printers
    void printVariationalParameters(int i);
};

#endif // TWOELECTRONJASTROV_H
