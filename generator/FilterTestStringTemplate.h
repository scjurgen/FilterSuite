
TEST(FilterCharacter, __PASSMODE____CHARACTER____ORDER__)
{
    for (double a = 0.0001; a < 0.5; a *= 10)
    {
        FilterCoefficients fc(__PASSMODE__, __CHARACTER__, __ORDER__);
        fc.chebyshevRipple = -3;
        fc.qfactor         = 10;
        fc.computeFilter(a, 1.3 * a);
        __FILTERNAME__<double> filter(__CHARACTER__);
        filter.setFilterParameters(fc);
        CalcBins(fc.getCurrentDisplaynameWithAlpha(), 100, &filter);
        CalcSample(fc.getCurrentDisplaynameWithAlpha(), &filter);
    }
}
