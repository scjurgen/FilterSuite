//
// Created by Jürgen Schwietering on 11/1/15.
//

#ifndef FILTERSUITE_FILTERCOEFFICIENTSPRINT_H
#define FILTERSUITE_FILTERCOEFFICIENTSPRINT_H

#include "FilterCoefficients.h"

class FilterCoefficientsPrint : public CalcFilter
    {
    public:
        FilterPrint (PASSMODE passMode = PM_BP, FILTERCHARACTER character=FC_BESSEL, int order=1)
                : CalcFilter(passMode,character,order)
        {

        }
        void printresults ()
        {
            printfilter();
            //if (opt_l)
            //{
            printf("G  = %.20e\n", hypot(pbgain.imag(), pbgain.real()));
            printcoeffs("NZ", zplane.zeros.size(), xcoeffs);
            printcoeffs("NP", zplane.poles.size(), ycoeffs);
        }

        void printcoeffs (const char *pz, int npz, const std::vector<double> &coeffs)
        {
            //echo "printcoeffs(pz, npz, coeffs)\n";
            printf("%s = %d\n", pz, npz);
            for (size_t i = 0; i <= npz; i++)
                printf("%18.20e\n", coeffs[i]);
        }

        void printfilter ()
        {

            //echo "printfilter()\n";
            printf("alpha low     = %14.20f\n", raw_alpha1);
            printf("alpha high    = %14.20f\n", raw_alpha2);
/*        if (!((FC_RESONATOR == optCharacter) | opt_w | opt_z))
        {
            printf("warped alpha low  = %14.10f\n", warped_alpha1);
            printf("warped alpha high = %14.10f\n", warped_alpha2);
        }
        */
            printgain("dc    ", dc_gain);
            printgain("centre", fc_gain);
            printgain("hf    ", hf_gain);
            printgain("rGain ", rGain);
            puts("\n");
            if (!(FC_RESONATOR == optCharacter))
                printrat_s();
            printrat_z();
            printrecurrence();
        }

        void printgain (const char *str, const std::complex<double> &gain)
        {
            //echo "printgain(str, gain)\n";
            double r = hypot(gain.imag(), gain.real());
            printf("gain at %s:   mag = %15.20e", str, r);
            if (r > EPS)
                printf("   phase = %14.20f (atan2)/3.14159", atan2(gain.imag(), gain.real()) / M_PI);
            puts("\n");
        }

        void printrat_s ()  /* prS-plane poles and zeros */
        {
            printf("S-plane zeros:\n");
            printpz(splane.zeros, splane.zeros.size());
            printf("S-plane poles:\n");
            printpz(splane.poles, splane.poles.size());
        }

        void printrat_z ()  /* prZ-plane poles and zeros */
        {
            //echo "printrat_z()	/* prZ-plane poles and zeros */\n";
            printf("Z-plane zeros:\n");
            printpz(zplane.zeros, zplane.zeros.size());
            printf("Z-plane poles:\n");
            printpz(zplane.poles, zplane.poles.size());
        }

        void printpz (const std::vector<std::complex<double> > &pzvec, int num)
        {
            //echo "printpz(pzvec, num)\n";
            int n1 = 0;
            while (n1 < num)
            {
                prcomplex(pzvec[n1]);
                int n2 = n1 + 1;
                while (n2 < num && pzvec[n2] == pzvec[n1])
                    n2++;
                if (n2 - n1 > 1)
                    printf("\t%d times", n2 - n1);
                n1 = n2;
            }
            printf("\n");
        }

        void printrecurrence ()
        {
            printf("Recurrence relation:\n");
            printf("y[n] = \n");


            for (size_t i = 0; i < zplane.zeros.size() + 1; i++) {
                double x = xcoeffs[i];
                if (x != 0) {
                    if (i > 0)
                        printf(" + ");
                    double f = fmod(abs(x), 1.0);
                    const char *fmt = (f < EPS || f > 1.0 - EPS) ? "%3g" : "%14.10f";
                    printf("(");
                    printf(fmt, x);
                    printf(" * x[n-%2d])\n", zplane.zeros.size() - i);
                }
            }
            puts("\n");
            for (size_t i = 0; i < zplane.poles.size(); i++)
            {
                printf("     + (%14.20f * y[n-%2d])\n", ycoeffs[i], zplane.poles.size() - i);
            }
            puts("\n");
        }

        void prcomplex (std::complex<double> z)
        {
            //echo "prcomplex(z)\n";
            printf("%14.20f + j %14.20f", z.real(), z.imag());
        }

    };
};


#endif //FILTERSUITE_FILTERCOEFFICIENTSPRINT_H
