//
// Created by juergen.schwietering@teufel.local on 11/6/15.
//

#ifndef FILTERSUITE_COMPLEXADDITIONS_H
#define FILTERSUITE_COMPLEXADDITIONS_H


#include <iostream>
#include <complex>
#include <vector>

class ComplexAdditions
{
public:
    inline std::complex<double> evaluatePolynomialInZ (
            const std::vector<std::complex<double> > &coeffs,
            std::complex<double> z)
    { /* evaluate polynomial in z, substituting for z */
        std::complex<double> sum = std::complex<double>(0.0, 0.0);
        for (int i = coeffs.size()-1; i >= 0; i--)
            sum = (sum * z) + coeffs[i];
        return sum;
    }

    inline std::complex<double> doubleComplexMultiply(double a, std::complex<double> z)
    {
        return std::complex<double>(a * z.real(), a * z.imag());
    }

    static double Xsqrt(double x)
    {
        /* because of deficiencies in hypot on Sparc, it's possible for arg of Xsqrt to be small and -ve,
         which logically it can't be (since r >= |x.re|).	 Take it as 0. */
        return (x >= 0.0) ? sqrt(x) : 0.0;
    }

protected:
    inline std::complex<double> bilinearTransform(std::complex<double> pz)
    {
        return  (2.0 + pz) / (2.0 - pz);
    }

    inline std::complex<double> cconj (std::complex<double> z)
    {
        return std::complex<double>(z.real(), -z.imag());
    }

    inline std::complex<double> complexExpj(double theta) // double
    {
        return std::complex<double>(cos(theta), sin(theta));
    }

    inline std::complex<double> complexcexp (std::complex<double> z)
    {
        return doubleComplexMultiply(exp(z.real()), complexExpj(z.imag()));
    }

    inline std::complex<double> csqrt(std::complex<double> x)
    {
        auto r = hypot(x.imag(), x.real());
        if (x.imag() >= 0.0)
        {
            return std::complex<double>(Xsqrt(0.5 * (r + x.real())),
                                        Xsqrt(0.5 * (r - x.real())));
        }
        return std::complex<double>(Xsqrt(0.5 * (r + x.real())),
                                    -Xsqrt(0.5 * (r - x.real())));
    }

    /**
             * evaluate response, substituting for z
             * */
    inline std::complex<double> evaluateResponse (
            const std::vector<std::complex<double> > &topco,
            const std::vector<std::complex<double> > &botco,
            std::complex<double> z)
    {

        std::complex<double> top = evaluatePolynomialInZ(topco, z);
        std::complex<double> bot = evaluatePolynomialInZ(botco, z);
        std::complex<double> res = top / bot;
        return res;
    }

    inline std::complex<double> reflect (std::complex<double> z)
    {
        //echo "reflect(z)\n";
        std::complex<double> r = std::complex<double>(hypot(z.imag(), z.real()), 0);
        return z / (r * r);
    }
};


#endif //FILTERSUITE_COMPLEXADDITIONS_H
