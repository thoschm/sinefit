

#include <SineFit.h>
#include <fstream>
#include <stdlib.h>


using namespace SineFit;


bool loadSequence(std::vector<float> *seq, const char *file)
{
    seq->clear();
    std::ifstream ifs;
    ifs.open(file, std::ios::in);
    if (ifs.fail())
    {
        std::cerr << "failed to read sequence!\n";
        return false;
    }
    float val;
    for ( ; ; )
    {
        ifs >> val;
        if (ifs.eof()) break;
        seq->push_back(val);
    }
    ifs.close();
    return true;
}


bool dumpSequence(const std::vector<float> &seq, const char *file)
{
    std::ofstream ofs;
    ofs.open(file, std::ios::out);
    if (ofs.fail())
    {
        std::cerr << "failed to write sequence!\n";
        return false;
    }
    for (uint i = 0; i < seq.size(); ++i)
    {
        ofs << seq[i] << std::endl;
    }
    ofs.close();
    return true;
}


#define WINDOW_SIZE 500


int main(int argc, char **argv)
{
    // check args
    if (argc != 2)
    {
        std::cerr << "Usage:\n   sinefit <sequence.txt>\n";
        return EXIT_FAILURE;
    }

    // load input sequence
    std::vector<float> indata, outdata;
    if (!loadSequence(&indata, argv[1]))
    {
        return EXIT_FAILURE;
    }
/*

    for (uint i = 0; i < 500u; ++i)
    {
        indata.push_back(std::sqrt(1 * i) + std::sin(0.1 * i) + std::sin(0.05 * (i + 17)) * std::cos(0.02 * (i + 23)) + 0.01f * i + 5.0f * std::sin(0.01f * (i + 100)));
    }
*/

    const uint limit = indata.size() - WINDOW_SIZE;
    outdata.resize(indata.size(), 0.0f);
    for (uint i = 0; i <= limit; i += 10u)
    {
        std::cerr << i << " / " << limit << std::endl;
        SineFitter<float> fitter(WINDOW_SIZE);
        SineParams<float> result = fitter.fit(indata, i, 0.0f, 2000);
        TurnResult turn = fitter.checkTurn(result, 0.9f);
        if (result.score < 0.1 && turn != NO_TURN) outdata[i + WINDOW_SIZE - 1u] = indata[i + WINDOW_SIZE - 1u];

        if (i == 1000u)
        {
            std::vector<float> gen;
            fitter.generateSine(&gen, indata, result, i);
            dumpSequence(gen, "gen.txt");
        }
    }

    dumpSequence(indata, "sine.txt");
    dumpSequence(outdata, "turn.txt");

/*
    SineFitter<float> fitter(1000);
    SineParams<float> result = fitter.fit(indata, UINT_MAX, 0.0f, 2000);

    std::vector<float> norm, gen;
    fitter.copyNormalized(&norm, indata, UINT_MAX);
    fitter.generateSine(&gen, indata, result, UINT_MAX);

    //std::cerr << result.freq << std::endl;

    dumpSequence(norm, "norm.txt");
    dumpSequence(gen, "gen.txt");

*/
    return 0;
}
