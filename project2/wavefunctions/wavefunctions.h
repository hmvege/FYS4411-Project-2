#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H


class WaveFunctions
{
protected:
    int nParticles;
    int nDimensions;
    int nVarParams;
    double coulomb(double **r);
public:
    WaveFunctions() { }
    WaveFunctions(int new_nParticles, int new_nDimensions, int new_nVarParams);
    virtual ~WaveFunctions() {}

    // Virtuals used by all other classes
    virtual double calculate(double **r);
    virtual double localEnergy(double **r);
    virtual void quantumForce(double **r, double **F, int k);
    virtual void steepestDescent(double **r, double E, double ESum, int NCycles);
    virtual void sampleSD(double **r, double E);
    // Printers
    virtual void printVariationalParameters();
};

#endif // WAVEFUNCTIONS_H
