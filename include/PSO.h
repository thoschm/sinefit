
#ifndef _PSO_H_
#define _PSO_H_


/////////////////////////////////
// INCLUDE
/////////////////////////////////
#include "XorShift.h"
#include <float.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <assert.h>


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
// Particle
/////////////////////////////////
template <typename NumericalType>
class Particle
{
public:
    Particle(const uint dim) : score((NumericalType)FLT_MAX),
                               tmp((NumericalType)FLT_MAX),
                               mDim(dim)
    {
        x.resize(dim, (NumericalType)FLT_MAX);
        v.resize(dim, (NumericalType)FLT_MAX);
        best.resize(dim, (NumericalType)FLT_MAX);
    }

    std::vector<NumericalType> x,
                               v,
                               best;
    NumericalType score,
                  tmp;

protected:
    uint mDim;
};


/////////////////////////////////
// PSO class
/////////////////////////////////
template <typename NumericalType>
class PSO
{
    PSO(const PSO &other);
    PSO &operator=(const PSO &other);

public:
    // alloc particles
    PSO(const uint dim,
        const uint particleCount,
        NumericalType (*scoreFunction)(const std::vector<NumericalType> *, const void *),
        const void *payload) :
        mDim(dim),
        mParticleCount(particleCount),
        scoreFunc(scoreFunction),
        mPayload(payload),
        mW((NumericalType)FLT_MAX),
        mCP((NumericalType)FLT_MAX),
        mCG((NumericalType)FLT_MAX),
        mBestScore((NumericalType)FLT_MAX)
    {
        // fix rnd
        mRnd.setSeed(0x12345678u);

        // alloc
        mBestPos.resize(dim, (NumericalType)FLT_MAX);
        mParticles.resize(particleCount, Particle<NumericalType>(dim));
    }

    // init
    void init(const std::vector<NumericalType> &lowerLimits,
              const std::vector<NumericalType> &upperLimits,
              const NumericalType w  = (NumericalType)0.7,
              const NumericalType cp = (NumericalType)1.4,
              const NumericalType cg = (NumericalType)1.4)
    {
        // check
        assert(lowerLimits.size() == mDim && upperLimits.size() == mDim);

        // params
        mW = w;
        mCP = cp;
        mCG = cg;

        // init particles
        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType> &par = mParticles[i];
            par.score = (NumericalType)FLT_MAX;
            par.tmp = (NumericalType)FLT_MAX;

            // init
            for (uint d = 0; d < mDim; ++d)
            {
                // x, v, best
                const NumericalType diff = upperLimits[d] - lowerLimits[d];
                par.x[d] = mRnd.uniform() * diff + lowerLimits[d];
                par.v[d] = mRnd.uniform() * (NumericalType)2.0 * diff - diff;
                par.best[d] = mRnd.uniform() * diff + lowerLimits[d];
            }
        }

        // reset swarm
        mBestScore = (NumericalType)FLT_MAX;
        std::fill(mBestPos.begin(), mBestPos.end(), (NumericalType)0.0);
    }

    // one optimization step
    NumericalType step()
    {
        // check
        assert(mW != (NumericalType)FLT_MAX);

        // compute scores
        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType> &par = mParticles[i];

            // compute cost
            par.tmp = scoreFunc(&(par.x), mPayload);
        }

        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType> &par = mParticles[i];

            // update scores
            if (par.tmp < par.score)
            {
                par.score = par.tmp;
                par.best = par.x;

                // swarm
                if (par.tmp < mBestScore)
                {
                    mBestScore = par.tmp;
                    mBestPos = par.x;
                }
            }

            // update position x
            for (uint d = 0; d < mDim; ++d)
            {
                // uniform random values
                const NumericalType rp = mRnd.uniform(),
                                    rg = mRnd.uniform();
                // update velocity and pos
                par.v[d] = mW * par.v[d] +
                           mCP * rp * (par.best[d] - par.x[d]) +
                           mCG * rg * (mBestPos[d] - par.x[d]);
                par.x[d] += par.v[d];
            }
        }
        return mBestScore;
    }

    // get particle pointer
    const Particle<NumericalType> &getParticles()
    {
        return mParticles;
    }

    // get current best
    const std::vector<NumericalType> &getBest()
    {
        return mBestPos;
    }

    NumericalType getScore()
    {
        return mBestScore;
    }

    // print particles
    void print()
    {
        std::cerr << "best score: " << mBestScore << std::endl;
        std::cerr << "best pos..: ";
        for (uint i = 0; i < mDim; ++i)
        {
            std::cerr << mBestPos[i] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "particles.: " << std::endl;
        for (uint i = 0; i < mParticleCount; ++i)
        {
            const Particle<NumericalType> &par = mParticles[i];
            for (uint d = 0; d < mDim; ++d)
            {
                std::cerr << par.x[d] << " ";
            }
            std::cerr << std::endl;
        }
        std::cerr << std::endl;
    }

    // dump particles to file
    void dump(const char *file)
    {
        std::ofstream of;
        of.open(file, std::ios::out);
        for (uint i = 0; i < mParticleCount; ++i)
        {
            const Particle<NumericalType> &par = mParticles[i];
            for (uint d = 0; d < mDim; ++d)
            {
                of << par.x[d] << " ";
            }
            of << std::endl;
        }
        of.close();
    }

protected:
    uint mDim, mParticleCount;
    std::vector<Particle<NumericalType> > mParticles;
    NumericalType (*scoreFunc)(const std::vector<NumericalType> *, const void *);
    const void *mPayload;
    XorShift<NumericalType> mRnd;

    // Params
    NumericalType mW,  // velocity inertia weight
                  mCP, // personal best weight
                  mCG; // group best weight

    // Swarm
    std::vector<NumericalType> mBestPos;
    NumericalType mBestScore;
};

}

#endif
