

#include "../FilterCoefficients.h"
#include "../FilterSuite.h"
#include "FftBinCompare.h"
#include "CalcBins.h"


int main(void)
{
   {
      BandPassOrder1<double> filter(FC_CHEBYSHEV);
      filter.setChebyshevRipple(-3.0);
      CalcAlphaSweepSample("SweepBandpassButterworthOrder1", &filter);
   }
   {
      BandPassOrder1<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandpassButterworthOrder1", &filter);
   }
   {
      BandStopOrder1<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandstopButterworthOrder1", &filter);
   }
    {
        BandPassOrder1<double> filter(FC_BESSEL);
        CalcAlphaSweepSample("SweepBandpassBesselOrder1", &filter);
    }
    {
        BandPassResonatorOrder1<double> filter(FC_RESONATOR);
        filter.setQfactor(100);
        CalcAlphaSweepSample("SweepBandpassResonatorOrder1", &filter);
    }



   {
      BandPassOrder2<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandpassButterworthOrder2", &filter);
   }
   {
      BandStopOrder2<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandstopButterworthOrder2", &filter);
   }
   {
      BandPassOrder2<double> filter(FC_BESSEL);
      CalcAlphaSweepSample("SweepBandpassBesselOrder2", &filter);
   }
   {
      BandPassOrder3<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandpassButterworthOrder3", &filter);
   }
   {
      BandStopOrder3<double> filter(FC_BUTTERWORTH);
      CalcAlphaSweepSample("SweepBandstopButterworthOrder3", &filter);
   }
   {
      BandPassOrder3<double> filter(FC_BESSEL);
      CalcAlphaSweepSample("SweepBandpassBesselOrder3", &filter);
   }
}
