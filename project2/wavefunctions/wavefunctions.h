#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H


class WaveFunctions
{
protected:
    int nParticles;
    int nDimensions;
public:
    WaveFunctions() { }
    WaveFunctions(int new_nParticles, int new_nDimensions);
    virtual ~WaveFunctions() {}

    virtual double calculate(double **positions);
    virtual double localEnergy(double ** positions);
};

#endif // WAVEFUNCTIONS_H
