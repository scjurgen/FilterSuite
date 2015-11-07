/*   This program is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; version 2 of the License.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 */

/* Original Code written by A.J. Fisher as defined below.  New
 * and updates peformed by mkfilter project team as hosted and
 * defined on http://mkfilter.sourceforge.net/
 */

/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Header file */
#ifndef MKFILTER_H
#define MKFILTER_H

#include <cstring>
#include <cstdlib>

#define global
#define unless(x)   if(!(x))
#define until(x)    while(!(x))

#define VERSION	    "4.6"

#ifndef	PI
    #define PI	    3.14159265358979323846264338327950288  /* Microsoft C++ does not define M_PI ! */
#endif

#define TWOPI	    (2.0 * PI)
#define EPS	    1e-10
#define MAXORDER    10
#define MAXPZ	    512	    /* .ge. 2*MAXORDER, to allow for doubling of poles in BP filter;
                               high values needed for FIR filters */
#define MAXSTRING   256

typedef void (*proc)();
typedef unsigned int uint;

extern "C"
{ 
    double atof(const char*);
    void exit(int);
};

extern char *progname;
extern void readdata(char*, double&, int&, double*, int&, double*);

inline double sqr(double x)	    { return x*x;			       }
inline bool seq(char *s1, char *s2) { return strcmp(s1,s2) == 0;	       }
inline bool onebit(uint m)	    { return (m != 0) && ((m & m-1) == 0);     }

inline double asinh(double x)
{ 
    /* Microsoft C++ does not define */
    return log(x + sqrt(1.0 + sqr(x)));
}

inline double fix(double x)
{ 
    /* nearest integer */
    return (x >= 0.0) ? floor(0.5+x) : -floor(0.5-x);
}

#endif /*MKFILTER_H*/
