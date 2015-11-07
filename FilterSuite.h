
#include <vector>

#include "FilterCoefficients.h"

template<typename _CalcType_>
class GenericFilter
{
public:
    virtual _CalcType_ step(_CalcType_) = 0;
    virtual void setAlpha(double alpha) = 0;
    virtual void setAlphas(double alphaLow, double alphaHigh) = 0;
};


#include "AllFilters.h"

