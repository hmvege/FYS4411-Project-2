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
public:
    twoElectronJastrov(int new_nParticles, int new_nDimensions, double new_omega, double new_alpha, double new_a, double new_beta);

    double calculate(double **r);
    double localEnergy(double **positions);
    void quantumForce(double **positions, double **F, int k);
    // Setters
    void setOmega(double newOmega) { omega = newOmega; }
    void setAlpha(double newAlpha) { alpha = newAlpha; }
    void setBeta(double newBeta) { beta = newBeta; }
    void set_a(double new_a) { a = new_a; }
};

#endif // TWOELECTRONJASTROV_H
