#ifndef METROPOLIS_H
#define METROPOLIS_H


class Metropolis
{
private:
    int nParticles;
    int nDimensions;

    double seed;

    void update();
    void sampleSystem();
    double R();

    // Temporary instance till more generalized class is made
    double (*waveFunction)(double ** rPos);
    double (*localEnergy)(double ** rPos);
public:
    Metropolis();
    ~Metropolis();
    void runMetropolis(unsigned int MCCycles);

    void setRNGSeed(double newSeed);

    // Setters
    void setNParticles(int new_nParticles) { nParticles = new_nParticles; }
    void setNDimensions(int new_nDimensions) { nDimensions= new_nDimensions; }

    // Temporary instance till more generalized class is made
    void setWaveFunction(double (*wf)(double ** rPos)) { waveFunction = wf; }
    void setEnergyFunc(double (*en)(double ** rPos)) { localEnergy= en;}
};

#endif // METROPOLIS_H
