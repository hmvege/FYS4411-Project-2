#include <random>
#include <iostream>
#include <ctime>
#include "metropolis.h"

using std::cout;
using std::endl;

Metropolis::Metropolis()
{

}

Metropolis::~Metropolis()
{

}

void Metropolis::update()
{
    double oldWaveFunction = 0;
    double newWaveFunction = 0;

    double ** oldPositions;
    double ** newPositions;

}

void Metropolis::sampleSystem()
{

}

double Metropolis::R()
{
    /*
     * Metropolis ratio
     */
}

void Metropolis::setRNGSeed(double newSeed)
{
    seed = newSeed;
}

void Metropolis::runMetropolis(unsigned int MCCycles)
{
    // Setting up random generators
    std::mt19937_64 generator(std::time(nullptr)); // Starting up the Mersenne-Twister19937 function
    std::uniform_real_distribution<double> uniform_distribution(0,1);
//    std::uniform_real_distribution<double> epsilon_distribution(-epsilon, epsilon);

    double stepLength = 1.0;

    double oldWaveFunction = 0;
    double newWaveFunction = 0;

    double E = 0;
    double ESum = 0;
    double ESumSquared = 0;

    double * MCSamples = new double[MCCycles];
    double ** rPositionsOld = new double * [nParticles];
    double ** rPositionsNew = new double * [nParticles];

    for (int i = 0; i < nParticles; i++)
    {
        rPositionsOld[i] = new double[nDimensions];
        rPositionsNew[i] = new double[nDimensions];
        for (int j = 0; j < nDimensions; j ++)
        {
            rPositionsOld[i][j] = stepLength * uniform_distribution(generator) - 0.5;
            rPositionsNew[i][j] = rPositionsOld[i][j];
        }
    }

    oldWaveFunction = waveFunction(rPositionsOld);
    for (unsigned int cycle = 0; cycle < MCCycles; cycle++)
    {
//        MCSamples[i] = 0;
        for (int i = 0; i < nParticles; i++)
        {
            for (int j = 0; j < nDimensions; j++)
            {
                rPositionsNew[i][j] = rPositionsOld[i][j] + stepLength * uniform_distribution(generator) - 0.5;
            }

            newWaveFunction = waveFunction(rPositionsNew); // Find the position with updated wavefunctions

            if (uniform_distribution(generator) <= (newWaveFunction*newWaveFunction)/(oldWaveFunction*oldWaveFunction))
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rPositionsOld[i][j] = rPositionsNew[i][j];
                    oldWaveFunction = newWaveFunction;
                }
            }
            else
            {
                for (int j = 0; j < nDimensions; j++)
                {
                    rPositionsNew[i][j] = rPositionsOld[i][j];
                }
            }
            E = localEnergy(rPositionsOld);
            ESum += E;
            ESumSquared += E*E;
        }
    }
    cout << "Energy = " << ESum/double(MCCycles * nParticles) << endl;
}
