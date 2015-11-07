/** This file is autogenerated, do not change manually
*/

#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-generated-matchers.h>
using ::testing::ElementsAre;
#include "../FilterCoefficients.h"
#include "../FilterSuite.h"
#include "FftBinCompare.h"
#include "CalcBins.h"


TEST(FilterCharacter, PM_LPFC_BESSEL1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BESSEL, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder1<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_BUTTERWORTH1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BUTTERWORTH, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder1<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_CHEBYSHEV1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_CHEBYSHEV, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder1<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BESSEL1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BESSEL, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder1<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BUTTERWORTH1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BUTTERWORTH, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder1<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_CHEBYSHEV1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_CHEBYSHEV, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder1<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BESSEL1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BESSEL, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder1<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BUTTERWORTH1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BUTTERWORTH, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder1<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_CHEBYSHEV1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_CHEBYSHEV, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder1<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BESSEL1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BESSEL, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder1<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BUTTERWORTH1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BUTTERWORTH, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder1<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_CHEBYSHEV1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_CHEBYSHEV, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder1<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_BESSEL2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BESSEL, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder2<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_BUTTERWORTH2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BUTTERWORTH, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder2<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_CHEBYSHEV2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_CHEBYSHEV, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder2<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BESSEL2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BESSEL, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder2<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BUTTERWORTH2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BUTTERWORTH, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder2<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_CHEBYSHEV2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_CHEBYSHEV, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder2<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BESSEL2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BESSEL, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder2<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BUTTERWORTH2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BUTTERWORTH, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder2<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_CHEBYSHEV2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_CHEBYSHEV, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder2<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BESSEL2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BESSEL, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder2<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BUTTERWORTH2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BUTTERWORTH, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder2<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_CHEBYSHEV2)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_CHEBYSHEV, 2);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder2<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_BESSEL3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BESSEL, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder3<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_BUTTERWORTH3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_BUTTERWORTH, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder3<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_LPFC_CHEBYSHEV3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_LP, FC_CHEBYSHEV, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        LowPassOrder3<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BESSEL3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BESSEL, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder3<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_BUTTERWORTH3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_BUTTERWORTH, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder3<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_HPFC_CHEBYSHEV3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_HP, FC_CHEBYSHEV, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        HighPassOrder3<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BESSEL3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BESSEL, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder3<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_BUTTERWORTH3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_BUTTERWORTH, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder3<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_CHEBYSHEV3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_CHEBYSHEV, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassOrder3<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BESSEL3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BESSEL, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder3<double> filter(FC_BESSEL);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_BUTTERWORTH3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_BUTTERWORTH, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder3<double> filter(FC_BUTTERWORTH);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_CHEBYSHEV3)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_CHEBYSHEV, 3);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopOrder3<double> filter(FC_CHEBYSHEV);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_APFC_RESONATOR1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_AP, FC_RESONATOR, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        AllPassResonatorOrder1<double> filter(FC_RESONATOR);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BPFC_RESONATOR1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BP, FC_RESONATOR, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandPassResonatorOrder1<double> filter(FC_RESONATOR);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}

TEST(FilterCharacter, PM_BSFC_RESONATOR1)
{
    for (double a=0.0001; a < 0.5; a*=10)
    {
        FilterCoefficients fc(PM_BS, FC_RESONATOR, 1);
        fc.chebyshevRipple = -3;
        fc.qfactor = 10;
        fc.computeFilter(a, 1.3*a);
        BandStopResonatorOrder1<double> filter(FC_RESONATOR);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}
