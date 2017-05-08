#ifndef IMPORTANCESAMPLER_H
#define IMPORTANCESAMPLER_H

#include "metropolisratio.h"
#include "../wavefunctions/wavefunctions.h"

class ImportanceSampler : public MetropolisRatio
{
private:
    // SHOULD THIS CLASS BE MERGED IWTH THE WAVEFUNCTIONS-INSTANCE??
    int NPart;
    int NDim;
    double D;
    double deltat;
    double deltatD;
    double exp_denom_factor;
    double denom_factor;

    WaveFunctions *WF = nullptr; // Will this create a double instance of the wavefunction, as we have one stored in the vmc?
public:
    ImportanceSampler();
    ImportanceSampler(double newD, double newDeltat, int newNPart, int newNDim);

    double G(double **rPosNew, double **rPosOld);
    double Ratio(double **rPosNew, double **rPosOld, double newWF, double oldWF);
    void setWaveFunction(WaveFunctions *newWF) { WF = newWF; }
};

#endif // IMPORTANCESAMPLER_H
