/*   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 of the License.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 */

/* Original Code was written by A.J. Fisher as defined below.
 * this adaptation written by Jürgen Schwietering (Berlin)
 */

/* gencode - use "-l" output from mkfilter to generate C code
   to implement digital filter
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   April 1995
*/

#ifndef FILTERSUITE_FILTERCOEFFICIENTS_H
#define FILTERSUITE_FILTERCOEFFICIENTS_H

#include <iostream>
#include <complex>
#include <vector>
#include <assert.h>

#include "ComplexAdditions.h"

typedef enum{
    PM_LP,
    PM_HP,
    PM_BP,
    PM_BS,
    PM_AP //allpass
} PASSMODE;

typedef enum
{
    FC_BESSEL,
    FC_BUTTERWORTH,
    FC_CHEBYSHEV,
    FC_RESONATOR,
    FC_PROPORTIONAL_INTEGRAL
}FILTERCHARACTER;


using namespace std;
//

struct PolesAndZerosRepresentation
{
public:
    std::vector< std::complex<double> > poles;
    std::vector< std::complex<double> > zeros;

    PolesAndZerosRepresentation()
    {
        clear();
    }

    void clear()
    {
        poles.clear();
        zeros.clear();
    }
};


class FilterCoefficients
    : public ComplexAdditions
{
public:
    const double EPS;

private:
    double rawAlphaLow;
    double rawAlphaHigh;
    unsigned int polemask;
    FILTERCHARACTER optCharacter;
    PASSMODE optPassmode;

private:
    bool usingPreWarp;  // -w		don't pre-warp
    bool usingMatchedZTransform;  // -z		use matched z-transform

    bool infq;
    size_t order;

    double warpedAlphaLow;
    double warpedAlphaHigh;

    PolesAndZerosRepresentation splane;
    PolesAndZerosRepresentation zplane;
    std::complex<double> dc_gain;
    std::complex<double> fc_gain;
    std::complex<double> hf_gain;
    std::complex<double> pbgain;
    std::vector<std::complex<double> > topcoeffs;
    std::vector<std::complex<double> > botcoeffs;

public:
    double rGain;
    double chebyshevRipple;
    double qfactor;
    std::vector<double> xcoeffs;
    std::vector<double> ycoeffs;

public:
    FilterCoefficients (PASSMODE passMode = PM_BP, FILTERCHARACTER character=FC_BESSEL, size_t order=1, double extra=0.0);

    void clearFilter();


public:
    bool isUsingPreWarp() const;

    void setUsingPreWarp(bool usingPreWarp);

    bool isUsingMatchedZTransform() const;

    void setUsingMatchedZTransform(bool usingMatchedZTransform);

    void computeFilter(double alow);

    void computeFilter(double alow, double ahigh);

    /* these should disappear from here
     */
    std::string getCurrentDisplayname ();

    static std::string getDisplayname (PASSMODE pm, FILTERCHARACTER filtercharacter, int order);
    std::string getCurrentDisplaynameWithAlpha ();
    std::string getPureFilterName(PASSMODE pm, FILTERCHARACTER filtercharacter, int order);
    std::string getCurrentPureFilterName();

private:
    void computeSPlanePoles ();

    void choosePole (std::complex<double> z);

    void prewarp ();

    void normalize ();

    void computeZByBilinearTransform ();

    void computeZPlaneByMatchedZTransform ();

    void resonatorComputeNotch ();

    void resonatorComputeAllPass ();

    void resonatorComputeBandPass ();

    void expandpoly ();

    void productOfPointsAsPolynomialOfZ (std::vector<std::complex<double> > &pz,
                                         std::vector<std::complex<double> > &coeffs);

    void multiplyFactorIntoCoefficents (std::complex<double> w, std::vector<std::complex<double> > &coeffs);

    void getGain ();

};

#endif //FILTERSUITE_FILTERCOEFFICIENTS_H
