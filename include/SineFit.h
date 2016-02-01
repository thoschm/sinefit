
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
    SineFitter(const uint windowSize,
               const uint particles,
               const NumericalType W,
               const NumericalType cp,
               const NumericalType cg) : mWindow(windowSize),
                                         mParticles(particles),
                                         mW(W),
                                         mCp(cp),
                                         mCg(cg)
    { }

    // compute min max normalization
    void normalize(const std::vector<NumericalType> &data,
                   const uint startAt,
                   NumericalType *minVal,
                   NumericalType *scaling)
    {
        // normalize data window to 0-1
        NumericalType vmin = (NumericalType)FLT_MAX,
                      vmax = (NumericalType)-FLT_MAX;
        for (uint k = 0; k < mWindow; ++k)
        {
            const uint idx = startAt + k;
            if (data[idx] < vmin) vmin = data[idx];
            if (data[idx] > vmax) vmax = data[idx];
        }
        NumericalType scale;
        if (vmin == vmax)
        {
            scale = (NumericalType)1.0;
            vmin -= (NumericalType)0.5;
        }
        else
        {
            scale = (NumericalType)1.0 / (vmax - vmin);
        }

        // commit
        *minVal = vmin;
        *scaling = scale;
    }

    // actual pso fitting
    SineParams<NumericalType> fit(const std::vector<NumericalType> &data)
    {
        SineParams<NumericalType> result;






        return result;
    }

private:
    uint mWindow,
         mParticles;
    NumericalType mW,
                  mCp,
                  mCg;
};

}

#endif
