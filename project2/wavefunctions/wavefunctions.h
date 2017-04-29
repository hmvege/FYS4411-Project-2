#ifndef WAVEFUNCTIONS_H
#define WAVEFUNCTIONS_H


class WaveFunctions
{
public:
    WaveFunctions();
    ~WaveFunctions();

    virtual calculateWF(double ** positions);
};

#endif // WAVEFUNCTIONS_H
