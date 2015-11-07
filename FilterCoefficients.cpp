//
// Created by Jürgen Schwietering on 11/1/15.
//

#include "FilterCoefficients.h"

static const std::vector< std::complex<double> > besselPoles =
    {
        {-1.00000000000e+00, 0.00000000000e+00},
        {-1.10160133059e+00, 6.36009824757e-01},
        {-1.32267579991e+00, 0.00000000000e+00},
        {-1.04740916101e+00, 9.99264436281e-01},
        {-1.37006783055e+00, 4.10249717494e-01},
        {-9.95208764350e-01, 1.25710573945e+00},
        {-1.50231627145e+00, 0.00000000000e+00},
        {-1.38087732586e+00, 7.17909587627e-01},
        {-9.57676548563e-01, 1.47112432073e+00},
        {-1.57149040362e+00, 3.20896374221e-01},
        {-1.38185809760e+00, 9.71471890712e-01},
        {-9.30656522947e-01, 1.66186326894e+00},
        {-1.68436817927e+00, 0.00000000000e+00},
        {-1.61203876622e+00, 5.89244506931e-01},
        {-1.37890321680e+00, 1.19156677780e+00},
        {-9.09867780623e-01, 1.83645135304e+00},
        {-1.75740840040e+00, 2.72867575103e-01},
        {-1.63693941813e+00, 8.22795625139e-01},
        {-1.37384121764e+00, 1.38835657588e+00},
        {-8.92869718847e-01, 1.99832584364e+00},
        {-1.85660050123e+00, 0.00000000000e+00},
        {-1.80717053496e+00, 5.12383730575e-01},
        {-1.65239648458e+00, 1.03138956698e+00},
        {-1.36758830979e+00, 1.56773371224e+00},
        {-8.78399276161e-01, 2.14980052431e+00},
        {-1.92761969145e+00, 2.41623471082e-01},
        {-1.84219624443e+00, 7.27257597722e-01},
        {-1.66181024140e+00, 1.22110021857e+00},
        {-1.36069227838e+00, 1.73350574267e+00},
        {-8.65756901707e-01, 2.29260483098e+00}
    };


FilterCoefficients::FilterCoefficients (PASSMODE passMode, FILTERCHARACTER character, size_t order, double extra)
        : EPS(1.0e-8)
        , optCharacter(character)
        , optPassmode(passMode)
        , order(order) {
    if (character==FC_CHEBYSHEV)
    {
        chebyshevRipple = extra;
    }
    if (character==FC_RESONATOR)
    {
        qfactor = extra;
    }
    usingPreWarp = true;
    usingMatchedZTransform = false;
    infq = false;
    rawAlphaLow = 0.01;
    rawAlphaHigh = 0.02;
    clearFilter();
}

void FilterCoefficients::clearFilter() {
    xcoeffs.clear();
    ycoeffs.clear();
    topcoeffs.clear();
    botcoeffs.clear();
    splane.clear();
    zplane.clear();
}

bool FilterCoefficients::isUsingPreWarp() const {
    return usingPreWarp;
}

void FilterCoefficients::setUsingPreWarp(bool usingPreWarp) {
    FilterCoefficients::usingPreWarp = usingPreWarp;
}

bool FilterCoefficients::isUsingMatchedZTransform() const {
    return usingMatchedZTransform;
}

void FilterCoefficients::setUsingMatchedZTransform(bool usingMatchedZTransform) {
    FilterCoefficients::usingMatchedZTransform = usingMatchedZTransform;
}

void FilterCoefficients::computeFilter(double alow) {
    computeFilter(alow,alow);
}

void FilterCoefficients::computeFilter(double alow, double ahigh) {
    rawAlphaLow = alow;
    rawAlphaHigh = ahigh;
    polemask = 0x7fffff;
    if (FC_RESONATOR == optCharacter)
    {
        assert(qfactor > 0.0);
        assert(qfactor < 100000.0);
        if (PM_BP == optPassmode) resonatorComputeBandPass();     /* bandpass resonator	 */
        else if (PM_BS == optPassmode) resonatorComputeNotch();     /* bandstop resonator (notch) */
        else if (PM_AP == optPassmode) resonatorComputeAllPass();     /* allpass resonator		 */
        else
        {
            printf("resonator is not defined for low or highpass");
            return;
        }
    }
    else
    {
        if (FC_PROPORTIONAL_INTEGRAL == optCharacter)
        {
            prewarp();
            splane.poles.push_back(std::complex<double>(0.0, 0.0));
            splane.zeros.push_back(std::complex<double>(-2.0 * M_PI * warpedAlphaLow, 0));
        }
        else
        {
            computeSPlanePoles();
            prewarp();
            normalize();
        }
        if (usingMatchedZTransform)
            computeZPlaneByMatchedZTransform();
        else
            computeZByBilinearTransform();
    }
    expandpoly();
    getGain();
}

void FilterCoefficients::computeSPlanePoles ()
{
    splane.clear();
    if (FC_BESSEL == optCharacter)
    { /* Bessel filter */
        //echo "bessel\n";
        size_t p = (order * order) / 4; /* ptr into table */
        if (order & 1)
        {
            choosePole(besselPoles[p]);
            p++;
        }
        for (size_t i = 0; i < order / 2; i++)
        {
            choosePole(besselPoles[p]);
            choosePole(cconj(besselPoles[p]));
            p++;
        }
    }
    if ((FC_BUTTERWORTH == optCharacter) | (FC_CHEBYSHEV == optCharacter))
    {
        for (size_t i = 0; i < 2 * order; i++)
        {
            auto theta = (order & 1) ? (i * M_PI) / order : ((i + 0.5) * M_PI) / order;
            choosePole(complexExpj(theta));
        }
    }
    if (FC_CHEBYSHEV == optCharacter)
    {
        // modify for Chebyshev (p. 136 DeFatta et al.)
        if (chebyshevRipple >= 0.0)
        {
            throw("Chebyshev ripple is positiv but must be less than 0.0\n");
        }
        auto rip = pow(10.0, -chebyshevRipple / 10.0);
        auto epsilon = sqrt(rip - 1.0);
        auto y = asinh(1.0 / epsilon) / (double)order;
        if (y <= 0.0)
        {
            printf("mkfilter: bug: Chebyshev y=%g; must be .gt. 0.0\n", y);
            exit(1);
        }
        for (size_t i = 0; i < splane.poles.size(); i++)
        {
            splane.poles[i] = std::complex<double>(splane.poles[i].real() * sinh(y),
                                                   splane.poles[i].imag() * cosh(y));
        }
    }
}

void FilterCoefficients::choosePole (std::complex<double> z)
{
    if (z.real() < 0.0)
    {
        if (polemask & 1)
            splane.poles.push_back(z);
        polemask >>= 1;
    }
}

void FilterCoefficients::prewarp ()
{

    if ((!usingPreWarp | usingMatchedZTransform))
    {
        warpedAlphaLow = rawAlphaLow;
        warpedAlphaHigh = rawAlphaHigh;
    }
    else
    {
        warpedAlphaLow = tan(M_PI * rawAlphaLow) / M_PI;
        warpedAlphaHigh = tan(M_PI * rawAlphaHigh) / M_PI;
    }
}

void FilterCoefficients::normalize ()
{
    auto w1 = 2.0 * M_PI * warpedAlphaLow;
    auto w2 = 2.0 * M_PI * warpedAlphaHigh;
    switch(optPassmode) {
        case PM_LP:
            for (size_t i = 0; i < splane.poles.size(); i++)
                splane.poles[i] = splane.poles[i] * w1;
            splane.zeros.clear();
            break;
        case PM_HP:
            for (size_t i = 0; i < splane.poles.size(); i++) {
                splane.poles[i] = w1 / splane.poles[i];
            }
            splane.zeros.clear();
            for (size_t i = 0; i < splane.poles.size(); i++) {
                splane.zeros.push_back(0.0);   /* also N zeros at (0,0) */
            }
            break;
        case PM_BP: {
            double w0 = sqrt(w1 * w2);
            double bw = w2 - w1;
            splane.poles.resize(splane.poles.size() * 2);
            splane.zeros.resize(splane.poles.size() / 2);
            for (size_t i = 0; i < splane.poles.size() / 2; i++) {
                auto hba = 0.5 * (splane.poles[i] * bw);
                auto t2 = w0 / hba;
                auto temp = csqrt(1.0 - t2 * t2);
                splane.poles[i] = hba * (1.0 + temp);
                splane.poles[i + splane.poles.size() / 2] = hba * (1.0 - temp);
            }
            for (size_t i = 0; i < splane.poles.size() / 2; i++)
                splane.zeros[i] = 0.0;   // also N zeros at (0,0)
        }
            break;
        case PM_BS: {
            auto w0 = sqrt(w1 * w2);
            auto bw = w2 - w1;
            splane.poles.resize(splane.poles.size() * 2);
            splane.zeros.resize(splane.poles.size());
            for (size_t i = 0; i < splane.poles.size() / 2; i++) {
                auto hba = 0.5 * bw / splane.poles[i];
                std::complex<double> t2 = std::complex<double>(w0, 0) / hba;
                std::complex<double> temp = std::sqrt(1.0 - t2 * t2);
                splane.poles[i] = hba * (1.0 + temp);
                splane.poles[i + splane.poles.size() / 2] = hba * (1.0 - temp);
            }
            for (size_t i = 0; i < splane.poles.size() / 2; i++)     // also 2N zeros at (0, +-w0)
            {
                splane.zeros[i] = std::complex<double>(0.0, w0);
                splane.zeros[i + splane.poles.size() / 2] = std::complex<double>(0.0, -w0);
            }
        }
            break;
    }
}

void FilterCoefficients::computeZByBilinearTransform ()
{
        //echo "computeZByBilinearTransform() /* given S-plane poles & zeros, compute Z-plane poles & zeros, by bilinear transform */\n";
        zplane.poles.resize(splane.poles.size());
        zplane.zeros.resize(splane.zeros.size());

        for (size_t i = 0; i < zplane.poles.size(); i++)
        {
            zplane.poles[i] = bilinearTransform(splane.poles[i]);
        }

        for (size_t i = 0; i < zplane.zeros.size(); i++)
            zplane.zeros[i] = bilinearTransform(splane.zeros[i]);
        while (zplane.zeros.size() < zplane.poles.size())
            zplane.zeros.push_back(-1.0);
    }

void FilterCoefficients::computeZPlaneByMatchedZTransform ()
{

    zplane.poles.resize(splane.poles.size());
    zplane.zeros.resize(splane.zeros.size());
    for (size_t i = 0; i < zplane.poles.size(); i++)
        zplane.poles[i] = complexcexp(splane.poles[i]);
    for (size_t i = 0; i < zplane.zeros.size(); i++)
        zplane.zeros[i] = complexcexp(splane.zeros[i]);
    //print_r(zplane);
}

void FilterCoefficients::resonatorComputeNotch ()
{

    resonatorComputeBandPass();    /* iterate to place poles */
    double theta = 2 * M_PI * rawAlphaLow;
    std::complex<double> zz = complexExpj(theta);  /* place zeros exactly */
    zplane.zeros[0] = zz;
    zplane.zeros[1] = cconj(zz);
}

void FilterCoefficients::resonatorComputeAllPass ()
{
    resonatorComputeBandPass();    // iterate to place poles
    zplane.zeros[0] = reflect(zplane.poles[0]);
    zplane.zeros[1] = reflect(zplane.poles[1]);
}

void FilterCoefficients::resonatorComputeBandPass ()
{
    zplane.poles.resize(2);
    zplane.zeros.resize(2);

    zplane.zeros[0] = 1.0;
    zplane.zeros[1] =-1.0;
    auto theta = 2 * M_PI * rawAlphaLow; // where we want the peak to be
    if (infq)
    { /* oscillator */
        std::complex<double> zp = complexExpj(theta);
        zplane.poles[0] = zp;
        zplane.poles[1] = cconj(zp);
    }
    else
    { /* must iterate to find exact pole positions */
        productOfPointsAsPolynomialOfZ(zplane.zeros, topcoeffs);
        auto r = exp(-theta / (2.0 * qfactor));
        auto thm = theta;
        auto th1 = 0.0;
        auto th2 = M_PI;
        auto cvg = false;
        for (size_t i = 0; i < 50 && !cvg; i++)
        {
            auto zp = r * complexExpj(thm);
            zplane.poles[0] = zp;
            zplane.poles[1] = cconj(zp);
            productOfPointsAsPolynomialOfZ(zplane.poles, botcoeffs);
            auto g = evaluateResponse(topcoeffs, botcoeffs, complexExpj(theta));
            auto phi = g.imag() / g.real(); /* approx to atan2 */
            if (phi > 0.0)
                th2 = thm;
            else
                th1 = thm;
            if (abs(phi) < EPS)
                cvg = true;
            thm = 0.5 * (th1 + th2);
        }
        if (!cvg)
            printf("<font color=\"red\">warning: failed to converge</font>\n");
        //print_r(zplane);
        //print_r(topcoeffs);
    }
}

void FilterCoefficients::expandpoly ()
{
    productOfPointsAsPolynomialOfZ(zplane.zeros, topcoeffs);
    productOfPointsAsPolynomialOfZ(zplane.poles, botcoeffs);

    dc_gain = evaluateResponse(topcoeffs, botcoeffs, 1.0);
    double theta = 2 * M_PI * 0.5 * (rawAlphaLow + rawAlphaHigh); /* "jwT" for centre freq. */
    fc_gain = evaluateResponse(topcoeffs, botcoeffs, complexExpj(theta));
    hf_gain = evaluateResponse(topcoeffs, botcoeffs, -1.0);
    xcoeffs.resize(zplane.zeros.size()+1);
    for (size_t i = 0; i <= zplane.zeros.size(); i++) {
        xcoeffs[i] = +(topcoeffs[i].real() / botcoeffs[zplane.poles.size()].real());
    }
    ycoeffs.resize(zplane.poles.size()+1);
    for (size_t i = 0; i <= zplane.poles.size(); i++) {
        ycoeffs[i] = -(botcoeffs[i].real() / botcoeffs[zplane.poles.size()].real());
    }
}

void FilterCoefficients::productOfPointsAsPolynomialOfZ (std::vector<std::complex<double> > &pz,
                                                         std::vector<std::complex<double> > &coeffs)
{
    coeffs.resize(pz.size()+1);
    coeffs[0] = 1.0;
    for (size_t i = 0; i < pz.size(); i++)
        coeffs[i + 1] = 0.0;
    for (size_t i = 0; i < pz.size(); i++)
        multiplyFactorIntoCoefficents(pz[i], coeffs);
    // check computed coeffsX of z^k are all real
    for (size_t i = 0; i < pz.size() + 1; i++)
    {
        if (abs(coeffs[i].imag()) > EPS)
        {
            throw("coeff of z is not real; poles/zeros are not complex conjugates\n");
        }
    }
}

void FilterCoefficients::multiplyFactorIntoCoefficents (std::complex<double> w, std::vector<std::complex<double> > &coeffs)
{
    auto nw = std::complex<double>(-w.real(), -w.imag());
    for (int i = coeffs.size()-1; i >= 1; i--)
        coeffs[i] = nw * coeffs[i] + coeffs[i - 1];
    coeffs[0] = nw * coeffs[0];
}

void FilterCoefficients::getGain ()
{
    pbgain = std::complex<double>(1.0, 0.0);
    if (FC_PROPORTIONAL_INTEGRAL == optCharacter) pbgain = hf_gain;
    else if (PM_LP == optPassmode) pbgain = dc_gain;
    else if (PM_HP == optPassmode) pbgain = hf_gain;
    else if ((PM_BP == optPassmode) | (PM_AP == optPassmode)) pbgain = fc_gain;
    else if (PM_BS == optPassmode) pbgain = std::sqrt(dc_gain * hf_gain);
    rGain = hypot(pbgain.imag(), pbgain.real());
}

std::string FilterCoefficients::getDisplayname (PASSMODE pm, FILTERCHARACTER filtercharacter, int order)
{

    std::stringstream name;

    switch(pm)
    {
        case PM_LP: name <<"LowPass"; break;
        case PM_HP: name <<"HighPass"; break;
        case PM_BP: name <<"BandPass"; break;
        case PM_BS: name <<"BandStop"; break;
        case PM_AP: name <<"AllPass"; break;
    }
    switch(filtercharacter)
    {
        case FC_BESSEL: name <<"Bessel"; break;
        case FC_BUTTERWORTH: name <<"Buttwerworth"; break;
        case FC_CHEBYSHEV: name <<"Chebyshev"; break;
        case FC_RESONATOR: name <<"Resonator"; break;
        case FC_PROPORTIONAL_INTEGRAL: name <<"ProportionalIntegral"; break;
    }
    name << "Order" << order;
    return name.str();
}

std::string FilterCoefficients::getPureFilterName (PASSMODE pm, FILTERCHARACTER filtercharacter, int order)
{
    std::stringstream name;

    switch(pm)
    {
        case PM_LP: name <<"LowPass"; break;
        case PM_HP: name <<"HighPass"; break;
        case PM_BP: name <<"BandPass"; break;
        case PM_BS: name <<"BandStop"; break;
        case PM_AP: name <<"AllPass"; break;
    }
    switch(filtercharacter)
    {
        case FC_BESSEL: name <<""; break;
        case FC_BUTTERWORTH: name <<""; break;
        case FC_CHEBYSHEV: name <<""; break;
        case FC_RESONATOR: name <<"Resonator"; break;
        case FC_PROPORTIONAL_INTEGRAL: name <<"ProportionalIntegral"; break;
    }
    name << "Order" << order;
    return name.str();
}

std::string FilterCoefficients::getCurrentDisplayname ()
{
    return getDisplayname(optPassmode, optCharacter, order);
}

std::string FilterCoefficients::getCurrentDisplaynameWithAlpha ()
{
    std::stringstream s;
    s << getDisplayname(optPassmode, optCharacter, order);
    s << "_" << rawAlphaLow;
    if (optPassmode == PM_BS || optPassmode == PM_BP)
        s << "_" << rawAlphaLow;
    return s.str();
}


std::string FilterCoefficients::getCurrentPureFilterName ()
{
    return getPureFilterName(optPassmode, optCharacter, order);
}