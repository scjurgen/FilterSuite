#pragma once

const double epsilon = 1E-8;

template<typename _CalcType_>
void AddTestSine(size_t samplesCount,
                 _CalcType_ samplerate,
                 _CalcType_ f,
                 _CalcType_ factor,
                 std::vector<_CalcType_> &result
)
{

    for (size_t i=0;i < samplesCount; ++i)
    {
        // calculate in double, then cast it
        _CalcType_ v = static_cast<_CalcType_>(factor*sin(static_cast<double>(i)*f/static_cast<double>(samplerate)*M_PI*2.0));
        result.push_back(v);
    }
}

template<typename _CalcType_>
void GenerateTestWhiteNoise(size_t samplesCount,
                            _CalcType_ factor,
                            std::vector<_CalcType_> &result
)
{
    result.resize(samplesCount);
    for (size_t i=0;i < samplesCount; ++i)
    {
        _CalcType_ v = -factor/2 + factor*static_cast<_CalcType_>(rand())/static_cast<_CalcType_>(RAND_MAX);
        result[i] = v;
    }
}

std::vector<double> CalcBins(const std::string &displayName, int runs, GenericFilter<double> *filter)
{

    FftBins fftBins(44100);
    for (int n=0; n < runs; ++n)
    {
        std::vector<double> wave;
        GenerateTestWhiteNoise<double>(16384, 1.0, wave);
        std::vector<double> waveResult;
        waveResult.resize(wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            waveResult[i] = filter->step(wave[i]);
        }
        FftBins::HammingWindow(waveResult);
        fftBins.AddSlice(waveResult);
    }

    std::stringstream pngFilename;
    pngFilename << displayName << ".png";
    FILE *fp;
    fp = fopen("data.gnuplot","w");
    if (fp)
    {
        fprintf(fp,"set terminal png nocrop enhanced size 1024,320 font \"arial,8\"\n");
        fprintf(fp,"set output '%s'\n", pngFilename.str().c_str());
        fprintf(fp,"set key bmargin center horizontal Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid\n");
        fprintf(fp,"set style data lines\n");
        fprintf(fp,"set title \"%s\"\n", displayName.c_str());
        fprintf(fp,"set xrange [ 0.00000 : 120.00000 ] noreverse nowriteback\n");
        fprintf(fp,"x = 0.0\n");
        fprintf(fp,"plot 'data.dat' using 1:2 title \"log\" with filledcurve\n");
        fprintf(fp,"\n");
        fclose(fp);
    }

    fp = fopen("data.dat","w");
    if (fp)
    {
        fprintf(fp, "# x y log(y)\n");
        auto bins = fftBins.getBins();
        fprintf(fp, "0.0 0 0\n");
        for (size_t i = 0; i < bins.size(); ++i)
        {
            fprintf(fp, "%f %g %g\n", (double) i+1, bins[i], log(bins[i]));
        }
        fprintf(fp, "%f 0 0\n",(double)bins.size()+1.0);
        fclose(fp);
    }

    system("gnuplot -c data.gnuplot");

    return fftBins.getBins(); // make a copy
}


void CalcSample(const std::string &displayName, GenericFilter<double> *filter)
{
        std::vector<double> wave;
        GenerateTestWhiteNoise<double>(163840, 1.0, wave);
        std::vector<double> waveResult;
        waveResult.resize(wave.size());
        for (size_t i = 0; i < wave.size(); ++i)
        {
            waveResult[i] = filter->step(wave[i]);
        }

    std::stringstream waveFilename;
    waveFilename << displayName << ".mp3";
    FILE *fp;
    fp = fopen("infile.raw","w");
    if (fp)
    {
        for (size_t i = 0; i < waveResult.size(); ++i)
        {
            int16_t v = waveResult[i]*32000.0;
            fwrite(&v,2,1,fp);
        }
        fclose(fp);
    }
    std::stringstream soxstr;
    soxstr << "sox -b 16 -e signed-integer -c 1 -r 44.1k -t raw infile.raw " << waveFilename.str();
    system(soxstr.str().c_str());
}


void CalcAlphaSweepSample(const std::string &displayName, GenericFilter<double> *filter)
{
    std::vector<double> wave;
    GenerateTestWhiteNoise<double>(44100*4, 1.0, wave);
    std::vector<double> waveResult;
    waveResult.resize(wave.size());

    size_t stepSize = waveResult.size()/120;
    size_t nextSlizeUpdate = 0;
    int currentSlize = 100;
    for (size_t i = 0; i < wave.size(); ++i)
    {
        if (i==nextSlizeUpdate)
        {
            nextSlizeUpdate =  i+stepSize;
            double f = 440.0*pow(2.0,(currentSlize-69.0)/12.0); // frequency of A6 and basenote 69=A6
            --currentSlize;
            filter->setAlphas(f/44100.0,f*1.3/44100.0);
        }
        waveResult[i] = filter->step(wave[i]);
    }

    std::stringstream waveFilename;
    waveFilename << displayName << ".mp3";
    FILE *fp;
    fp = fopen("infile.raw","w");

    if (fp)
    {
        for (size_t i = 0; i < waveResult.size(); ++i)
        {
            int16_t v = waveResult[i]*32000.0;
            fwrite(&v,2,1,fp);
        }
        fclose(fp);
    }
    std::stringstream soxstr;
    soxstr << "sox -b 16 -e signed-integer -c 1 -r 44.1k -t raw infile.raw " << waveFilename.str();
    system(soxstr.str().c_str());
}

