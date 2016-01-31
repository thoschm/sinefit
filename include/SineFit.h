
#ifndef _SINE_FIT_H_
#define _SINE_FIT_H_


/////////////////////////////////
// INCLUDE
/////////////////////////////////
#include "PSO.h"


/////////////////////////////////
// TYPEDEFS
/////////////////////////////////
typedef unsigned int uint;


/////////////////////////////////
// NAMESPACE
/////////////////////////////////
namespace SineFit
{


/////////////////////////////////
// SINE FITTER
/////////////////////////////////
template <typename NumericalType>
struct SineParams
{
    NumericalType m, n, c,
                  amp, freq;
    bool success;
    SineParams() : m((NumericalType)0.0),
                   n((NumericalType)0.0),
                   c((NumericalType)0.0),
                   amp((NumericalType)0.0),
                   freq((NumericalType)0.0),
                   success(false)
    { }
};


template <typename NumericalType>
class SineFitter
{
private:
    SineFitter(const SineFitter &other);
    SineFitter &operator=(const SineFitter &other);

public:
    SineFitter(const uint particles,
               const NumericalType W,
               const NumericalType cp,
               const NumericalType cg) : mParticles(particles),
                                         mW(W),
                                         mCp(cp),
                                         mCg(cg)
    { }

    // actual pso fitting
    SineParams<NumericalType> fit(const std::vector<NumericalType> &data)
    {
        SineParams<NumericalType> result;






        return result;
    }

private:
    uint mParticles;
    NumericalType mW,
                  mCp,
                  mCg;
};

}

#endif
