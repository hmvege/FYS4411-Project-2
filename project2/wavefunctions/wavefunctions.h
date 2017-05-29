#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H


class WaveFunctions
{
protected:
    int nParticles;
    int nDimensions;
    int nVarParams;
    bool coulombInteraction = true;
    double coulomb(double **r);
    double r_ij(double *r1, double *r2);
    double SDStepLength = 0.01;
public:
    WaveFunctions() { }
    WaveFunctions(int new_nParticles, int new_nDimensions, int new_nVarParams);
    virtual ~WaveFunctions() { }

    // Virtuals used by all other classes
    virtual void initializeWFSampling(double **r);
    virtual double initializeWaveFunction(double **r);
//    virtual void initialize(double **r, double &WF);
    virtual double calculate(double **r, int k);
    virtual double localEnergy(double **r);
    virtual void quantumForce(double **r, double **F, int k);
    virtual void steepestDescent(double &E, int NCycles);
    virtual void sampleSD(double **r, double &E);
    virtual bool SDConvergenceCriteria();
    virtual void updateWF();
    virtual void revert();
    // Printers
    virtual void printVariationalParameters();
    // Setters
    void setSDStepLength(double newSDStepLength) { SDStepLength = newSDStepLength; }
    void setCoulombInteraction(bool C) { coulombInteraction = C; } // TEMP?
};

#endif // WAVEFUNCTIONS_H
