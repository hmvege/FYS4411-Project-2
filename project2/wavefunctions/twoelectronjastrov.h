#ifndef TWOELECTRONJASTROV_H
#define TWOELECTRONJASTROV_H

#include "wavefunctions.h"

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
    twoElectronJastrov(int new_nParticles, int new_nDimensions, int new_nVarParams, double new_omega, double new_alpha, double new_a, double new_beta);

//    void initialize(double **r, double &WF);
    void initializeWFSampling(double **r);
    double initializeWaveFunction(double **r);
    double calculate(double **r, int k);
    double localEnergy(double **r);
    void quantumForce(double **r, double **F, int k);
    void steepestDescent(double &ESum, int NCycles);
    void sampleSD(double **r, double &E);

    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }
    void set_a(double new_a) { a = new_a; }
    // Printers
    void printVariationalParameters();


//    // TEMP
//    void jastrov(bool t);
//    void coloumb(bool t);
};

#endif // TWOELECTRONJASTROV_H
