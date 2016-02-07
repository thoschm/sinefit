
#ifndef _SINE_FIT_H_
#define _SINE_FIT_H_


/////////////////////////////////
// INCLUDE
/////////////////////////////////
#include "PSO.h"
#include <limits.h>


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
                  amp, freq,
                  vmin, scale;
    bool success;
    SineParams() : m((NumericalType)0.0),
                   n((NumericalType)0.0),
                   c((NumericalType)0.0),
                   amp((NumericalType)0.0),
                   freq((NumericalType)0.0),
                   vmin((NumericalType)0.0),
                   scale((NumericalType)0.0),
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
               const uint particles = 100u,
               const NumericalType W =  (NumericalType)0.7,
               const NumericalType cp = (NumericalType)1.4,
               const NumericalType cg = (NumericalType)1.4) : mWindow(windowSize),
                                                              mParticles(particles),
                                                              mW(W),
                                                              mCp(cp),
                                                              mCg(cg)
    {
        const NumericalType wnd = windowSize;
        // prepare limits
        // m
        mLow.push_back(-1.0 / (wnd - (NumericalType)1.0));
        mHigh.push_back(1.0 / (wnd - (NumericalType)1.0));
        // n
        mLow.push_back((NumericalType)0.0);
        mHigh.push_back((NumericalType)1.0);
        // amp
        mLow.push_back((NumericalType)-1.0);
        mHigh.push_back((NumericalType)1.0);
        // freq
        const NumericalType halfperiod = M_PI / windowSize;
        mLow.push_back((NumericalType)1.0 * halfperiod);
        mHigh.push_back(wnd * halfperiod);
        // c
        mLow.push_back(-wnd);
        mHigh.push_back(wnd);
    }


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


    struct Payload
    {
        const std::vector<NumericalType> *sequence;
        uint elements;
    };


    // actual pso fitting
    SineParams<NumericalType> fit(const std::vector<NumericalType> &data,
                                  const uint startAt = UINT_MAX,
                                  const NumericalType breakScore = (NumericalType)0.0,
                                  const uint breakLoops = 1000u)
    {
        SineParams<NumericalType> result;
        const uint pos = std::min(startAt, (uint)(data.size() - mWindow));

        // get normalization parameters
        NumericalType minVal, scaling;
        normalize(data, pos, &minVal, &scaling);

        // copy window and normalize
        std::vector<NumericalType> window(mWindow);
        for (uint i = 0; i < mWindow; ++i)
        {
            window[i] = scaling * (data[pos + i] - minVal);
        }

        // fill payload
        Payload pl;
        pl.elements = mWindow;
        pl.sequence = &window;

        // init pso, TODO: move init into ctor
        PSO<NumericalType, 5> pso(mParticles, scoreFunc, (const void *)&pl);
        pso.init(mLow, mHigh);
        //pso.print();

        // optimize
        NumericalType s;
        uint l = 0;
        while ((s = std::sqrt(pso.step())) > breakScore)
        {
            std::cerr << "\roptimization error = " << s;
            if (++l == breakLoops) break;
            //pso.print();
        }
        std::cerr << "\roptimization error = " << s << std::endl;
        const NumericalType *values = pso.getBest();

        // fill return struct
        result.success = true;
        result.m = values[0];
        result.n = values[1];
        result.amp = values[2];
        result.freq = values[3];
        result.c = values[4];
        result.vmin = minVal;
        result.scale = scaling;
        return result;
    }


    // compute score for a single particle
    static NumericalType scoreFunc(const NumericalType *values, const uint dim, const void *payload)
    {
        // get payload
        const Payload *pl = (const Payload *)payload;
        // extract parameters
        const uint windowSize = pl->elements;
        const std::vector<NumericalType> *data = pl->sequence;
        const NumericalType m = values[0],
                            n = values[1],
                            amp = values[2],
                            freq = values[3],
                            c = values[4];
        // compute error
        NumericalType err = (NumericalType)0.0;
        for (uint i = 0; i < windowSize; ++i)
        {
            const NumericalType f = m * (NumericalType)i + n + amp * std::sin(freq * ((NumericalType)i + c)),
                                d = data->at(i) - f;
            err += d * d;
        }
        return err / windowSize;
    }

    // copy normalized data into vector
    void copyNormalized(std::vector<NumericalType> *out,
                        const std::vector<NumericalType> &data,
                        const uint startAt = UINT_MAX)
    {
        out->clear();
        out->resize(mWindow);

        // get normalization parameters
        NumericalType minVal, scaling;
        const uint pos = std::min(startAt, (uint)(data.size() - mWindow));
        normalize(data, pos, &minVal, &scaling);

        // copy window and normalize
        for (uint i = 0; i < mWindow; ++i)
        {
            out->at(i) = scaling * (data[pos + i] - minVal);
        }
    }

    // draw sine sequence
    void generateSine(std::vector<NumericalType> *out,
                      const SineParams<NumericalType> &p)
    {
        out->clear();
        out->resize(mWindow);

        // copy window and normalize
        for (uint i = 0; i < mWindow; ++i)
        {
            out->at(i) = p.m * (NumericalType)i + p.n + p.amp * std::sin(p.freq * ((NumericalType)i + p.c));
        }
    }

    // copy normalized data into vector
    void generateSine(std::vector<NumericalType> *out,
                      const std::vector<NumericalType> &data,
                      const SineParams<NumericalType> &p,
                      const uint startAt = UINT_MAX)
    {
        out->clear();
        out->resize(data.size(), (NumericalType)0.0);

        // get normalization parameters
        const NumericalType inv = (NumericalType)1.0 / p.scale;
        const uint pos = std::min(startAt, (uint)(data.size() - mWindow));

        // copy window and normalize
        for (uint i = 0; i < mWindow; ++i)
        {
            const NumericalType f = p.m * (NumericalType)i + p.n + p.amp * std::sin(p.freq * ((NumericalType)i + p.c));
            out->at(pos + i) = inv * f + p.vmin;
        }
    }

private:
    uint mWindow,
         mParticles;
    NumericalType mW,
                  mCp,
                  mCg;
    std::vector<NumericalType> mLow,
                               mHigh;
};

}

#endif
