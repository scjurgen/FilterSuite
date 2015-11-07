
#pragma once

#include <vector>
#include <fftw3.h>
#include <math.h>

#define FFTWIDTH 16384 

#define OCTAVEDIVISION 12

#define MAX_SPECTROGRAM_BINS 10*OCTAVEDIVISION


class FftBins
{
public:
    FftBins (size_t samplerate);
    ~FftBins ();

    void AddSlice(const std::vector<double> &data);
    size_t m_lastValidIndex;

    static void HammingWindow(std::vector<double> &data) {
        double omega = 2.0 * M_PI / (data.size()-1);
        for (size_t i=0; i < data.size(); ++i)
        {
            data[i] = data[i] * (0.54 - 0.46 * cos(omega* (i - 1)));
        }
    }
public:
    const std::vector<double> &getBins () const
    {
        return bins;
    }

private:
    size_t c_samplerate;
    std::vector<double> m_frequencySplits;
    fftw_plan m_fftplan;
    double m_in[FFTWIDTH];
    fftw_complex m_out[FFTWIDTH];
    std::vector<double> bins;

};


