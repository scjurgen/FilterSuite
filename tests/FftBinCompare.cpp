
#include "FftBinCompare.h"

#include <fftw3.h>
#include <math.h>

#include <sstream>
#include <iostream>


double maxval=0.0;

FftBins::FftBins (size_t samplerate)
    :  m_lastValidIndex(0)
    , c_samplerate(samplerate)
{
    m_frequencySplits.clear();


    m_fftplan = fftw_plan_dft_r2c_1d(FFTWIDTH, m_in, m_out, FFTW_MEASURE);

    // create a list of the frequency bands
    // there should be around 8-10 per octave, split by octavedivision
    // the last one will collect all frquencies until nyquist

    for (auto i=12.0; i < 128.0; i += 12.0/OCTAVEDIVISION)
    {
        double f = 440.0*pow(2,(i-69)/12); // frequency of A6 and basenote 69=A6
        m_frequencySplits.push_back(f);
    }
    m_frequencySplits.push_back(c_samplerate*2);
}

FftBins::~FftBins () {
    fftw_destroy_plan(m_fftplan);
}

void FftBins::AddSlice(const std::vector<double> &data)
{

    double power_spectrum[FFTWIDTH/2+1];

    for (size_t i = 0; i < FFTWIDTH;++i)
    {
        if (i >= data.size())
        {
            m_in[i] = 0;
        }
        else
        {
            m_in[i] = data[i];
        }
    }
    fftw_execute(m_fftplan);
    for (int k = 1; k < (FFTWIDTH + 1) / 2; ++k) {
        power_spectrum[k] = sqrt(m_out[k][0] * m_out[k][0] + m_out[k][1] * m_out[k][1]);
    }

    double fadvance = (double)c_samplerate/(double)FFTWIDTH;
    double f = fadvance;

    int currentFrequencyBin = 0;
    double ckcnt = 0;
    double sum = 0.0;
    bool newBins = bins.size()==0;
    for (int j=1; j < FFTWIDTH/2; ++j) {
        if (f>=m_frequencySplits[currentFrequencyBin]) {
            if (ckcnt>0.0) {
                if (newBins)
                    bins.push_back(sum / ckcnt);
                else
                    bins[currentFrequencyBin] += sum/ckcnt;

            }
            ckcnt = 0;
            sum = 0;
            currentFrequencyBin++;
        }
        sum += power_spectrum[j];
        ckcnt++;
        f += fadvance;
    }
    if (ckcnt>0.0) {
        if (newBins)
            bins.push_back(sum / ckcnt);
        else
            bins[currentFrequencyBin] += sum/ckcnt;
    }
}
