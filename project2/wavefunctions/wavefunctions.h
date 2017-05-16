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
    virtual double calculate(double **positions);
    virtual double localEnergy(double **positions);
    virtual void quantumForce(double **positions, double **F, int k);
    virtual void steepestDescent(double **rOld);
    // Printers
    virtual void printVariationalParameters();
};

#endif // WAVEFUNCTIONS_H
