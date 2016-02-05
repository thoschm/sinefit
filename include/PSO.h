
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
template <typename NumericalType, int Dim>
struct Particle
{
    NumericalType x[Dim],
                  v[Dim],
                  best[Dim],
                  score,
                  tmp;
};


/////////////////////////////////
// PSO class
/////////////////////////////////
template <typename NumericalType, int Dim>
class PSO
{
    PSO(const PSO &other);
    PSO &operator=(const PSO &other);

public:
    // alloc particles
    PSO(const uint particleCount,
        NumericalType (*scoreFunction)(const NumericalType *, const uint, const void *),
        const void *payload) :
        mParticleCount(particleCount),
        mParticles(NULL),
        scoreFunc(scoreFunction),
        mPayload(payload),
        mW((NumericalType)0.0),
        mCP((NumericalType)0.0),
        mCG((NumericalType)0.0),
        mBestScore((NumericalType)FLT_MAX)
    {
        // fix rnd
        //mRnd.setSeed(0x12345678u);

        memset(mBestPos, 0, Dim * sizeof(NumericalType));
        mParticles = new Particle<NumericalType, Dim>[particleCount];
    }

    // release particles
    ~PSO()
    {
        delete[] mParticles;
    }

    // init
    void init(const std::vector<NumericalType> &lowerLimits,
              const std::vector<NumericalType> &upperLimits,
              const NumericalType w  = (NumericalType)0.7,
              const NumericalType cp = (NumericalType)1.4,
              const NumericalType cg = (NumericalType)1.4)
    {
        // check
        if (lowerLimits.size() != Dim || upperLimits.size() != Dim)
        {
            std::cerr << "invalid limit vector!\n";
            return;
        }

        // params
        mW = w;
        mCP = cp;
        mCG = cg;

        // init particles
        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType, Dim> &par = mParticles[i];
            par.score = (NumericalType)FLT_MAX;
            par.tmp = (NumericalType)FLT_MAX;

            // init
            for (uint d = 0; d < Dim; ++d)
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
        memset(mBestPos, 0, Dim * sizeof(NumericalType));
    }

    // one optimization step
    NumericalType step()
    {
        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType, Dim> &par = mParticles[i];

            // compute cost
            par.tmp = scoreFunc(par.x, Dim, mPayload);
        }

        for (uint i = 0; i < mParticleCount; ++i)
        {
            // for each particle
            Particle<NumericalType, Dim> &par = mParticles[i];

            // update scores
            if (par.tmp < par.score)
            {
                par.score = par.tmp;
                memcpy(par.best, par.x, Dim * sizeof(NumericalType));

                // swarm
                if (par.tmp < mBestScore)
                {
                    mBestScore = par.tmp;
                    memcpy(mBestPos, par.x, Dim * sizeof(NumericalType));
                }
            }

            // update position x
            for (uint d = 0; d < Dim; ++d)
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
    const Particle<NumericalType, Dim> *getParticles()
    {
        return mParticles;
    }

    // get current best
    const NumericalType *getBest()
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
        for (uint i = 0; i < Dim; ++i)
        {
            std::cerr << mBestPos[i] << " ";
        }
        std::cerr << std::endl;
        std::cerr << "particles.: " << std::endl;
        for (uint i = 0; i < mParticleCount; ++i)
        {
            const Particle<NumericalType, Dim> &par = mParticles[i];
            for (uint d = 0; d < Dim; ++d)
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
            const Particle<NumericalType, Dim> &par = mParticles[i];
            for (uint d = 0; d < Dim; ++d)
            {
                of << par.x[d] << " ";
            }
            of << std::endl;
        }
        of.close();
    }

protected:
    uint mParticleCount;
    Particle<NumericalType, Dim> *mParticles;
    NumericalType (*scoreFunc)(const NumericalType *, const uint, const void *);
    const void *mPayload;
    XorShift<NumericalType> mRnd;

    // Params
    NumericalType mW,  // velocity inertia weight
                  mCP, // personal best weight
                  mCG; // group best weight

    // Swarm
    NumericalType mBestPos[Dim],
                  mBestScore;
};

}

#endif
