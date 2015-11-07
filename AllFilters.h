/** This file is autogenerated, do not change manually
*/

/**
 * Autogenerated filter LowPassOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class LowPassOrder1 : public GenericFilter<_CalcType_> {

public:
    LowPassOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_LP,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 2; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 1 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[1] = gain * x  + v[0] * ycoeffs[0];
        return v[0] + v[1];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[1];
    _CalcType_ ycoeffs[1];
    _CalcType_ v[2];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 1; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 1; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter HighPassOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class HighPassOrder1 : public GenericFilter<_CalcType_> {

public:
    HighPassOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_HP,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 2; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 1 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[1] = gain * x  + v[0] * ycoeffs[0];
        return -v[0] + v[1];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[1];
    _CalcType_ ycoeffs[1];
    _CalcType_ v[2];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 1; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 1; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandPassOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandPassOrder1 : public GenericFilter<_CalcType_> {

public:
    BandPassOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_BP,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return -v[0] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandStopOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandStopOrder1 : public GenericFilter<_CalcType_> {

public:
    BandStopOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_BS,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return v[0] * xcoeffs[0] + v[1] * xcoeffs[1] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter LowPassOrder2
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class LowPassOrder2 : public GenericFilter<_CalcType_> {

public:
    LowPassOrder2(FILTERCHARACTER fc)
    : filterCoefficients(PM_LP,fc,2) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return v[0] + 2*v[1] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter HighPassOrder2
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class HighPassOrder2 : public GenericFilter<_CalcType_> {

public:
    HighPassOrder2(FILTERCHARACTER fc)
    : filterCoefficients(PM_HP,fc,2) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return v[0] + -2*v[1] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandPassOrder2
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandPassOrder2 : public GenericFilter<_CalcType_> {

public:
    BandPassOrder2(FILTERCHARACTER fc)
    : filterCoefficients(PM_BP,fc,2) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 5; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 4 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[4] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2] + v[3] * ycoeffs[3];
        return v[0] + -2*v[2] + v[4];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[4];
    _CalcType_ ycoeffs[4];
    _CalcType_ v[5];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 4; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 4; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandStopOrder2
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandStopOrder2 : public GenericFilter<_CalcType_> {

public:
    BandStopOrder2(FILTERCHARACTER fc)
    : filterCoefficients(PM_BS,fc,2) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 5; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 4 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[4] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2] + v[3] * ycoeffs[3];
        return v[0] * xcoeffs[0] + v[1] * xcoeffs[1] + v[2] * xcoeffs[2] + v[3] * xcoeffs[3] + v[4];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[4];
    _CalcType_ ycoeffs[4];
    _CalcType_ v[5];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 4; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 4; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter LowPassOrder3
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class LowPassOrder3 : public GenericFilter<_CalcType_> {

public:
    LowPassOrder3(FILTERCHARACTER fc)
    : filterCoefficients(PM_LP,fc,3) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 4; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 3 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[3] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2];
        return v[0] + 3*v[1] + 3*v[2] + v[3];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[3];
    _CalcType_ ycoeffs[3];
    _CalcType_ v[4];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 3; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 3; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter HighPassOrder3
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class HighPassOrder3 : public GenericFilter<_CalcType_> {

public:
    HighPassOrder3(FILTERCHARACTER fc)
    : filterCoefficients(PM_HP,fc,3) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 4; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 3 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[3] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2];
        return -v[0] + 3*v[1] + -3*v[2] + v[3];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[3];
    _CalcType_ ycoeffs[3];
    _CalcType_ v[4];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 3; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 3; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandPassOrder3
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandPassOrder3 : public GenericFilter<_CalcType_> {

public:
    BandPassOrder3(FILTERCHARACTER fc)
    : filterCoefficients(PM_BP,fc,3) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 7; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 6 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[6] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2] + v[3] * ycoeffs[3] + v[4] * ycoeffs[4] + v[5] * ycoeffs[5];
        return -v[0] + 3*v[2] + -3*v[4] + v[6];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[6];
    _CalcType_ ycoeffs[6];
    _CalcType_ v[7];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 6; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 6; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandStopOrder3
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandStopOrder3 : public GenericFilter<_CalcType_> {

public:
    BandStopOrder3(FILTERCHARACTER fc)
    : filterCoefficients(PM_BS,fc,3) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 7; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 6 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[6] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1] + v[2] * ycoeffs[2] + v[3] * ycoeffs[3] + v[4] * ycoeffs[4] + v[5] * ycoeffs[5];
        return v[0] * xcoeffs[0] + v[1] * xcoeffs[1] + v[2] * xcoeffs[2] + v[3] * xcoeffs[3] + v[4] * xcoeffs[4] + v[5] * xcoeffs[5] + v[6];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[6];
    _CalcType_ ycoeffs[6];
    _CalcType_ v[7];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 6; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 6; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter AllPassResonatorOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class AllPassResonatorOrder1 : public GenericFilter<_CalcType_> {

public:
    AllPassResonatorOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_AP,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return v[0] * xcoeffs[0] + v[1] * xcoeffs[1] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandPassResonatorOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandPassResonatorOrder1 : public GenericFilter<_CalcType_> {

public:
    BandPassResonatorOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_BP,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return -v[0] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};


/**
 * Autogenerated filter BandStopResonatorOrder1
 * DO NOT CHANGE MANUALLY
 */
template<typename _CalcType_>
class BandStopResonatorOrder1 : public GenericFilter<_CalcType_> {

public:
    BandStopResonatorOrder1(FILTERCHARACTER fc)
    : filterCoefficients(PM_BS,fc,1) {
        reset();
    }

    void reset(void) {
        for (size_t i = 0; i < 3; ++i) {
            v[i] = 0;
        }
    }

    void setAlphas(double alphaLow, double alphaHigh) override {
        filterCoefficients.computeFilter(alphaLow, alphaHigh);
        setFilterParameters(1.0 / filterCoefficients.rGain, filterCoefficients.xcoeffs, filterCoefficients.ycoeffs);
    }

    void setAlpha(double alphaLow) override {
        setAlphas(alphaLow,alphaLow);
    }

    void setChebyshevRipple(double ripple) {
        filterCoefficients.chebyshevRipple = ripple;
    }

    void setQfactor(double q) {
        filterCoefficients.qfactor = q;
    }

    void inline shiftChain() {
        for (size_t j = 0; j < 2 ; ++j) {
            v[j] = v[j + 1];
        }
    }
    /** do not use in more than 1 call */
    _CalcType_ step(_CalcType_ x) override {
        shiftChain();
        v[2] = gain * x  + v[0] * ycoeffs[0] + v[1] * ycoeffs[1];
        return v[0] + v[1] * xcoeffs[1] + v[2];
    }

private:
    _CalcType_ gain;
    _CalcType_ xcoeffs[2];
    _CalcType_ ycoeffs[2];
    _CalcType_ v[3];
    FilterCoefficients filterCoefficients;

    void setFilterParameters(
            _CalcType_ gain,
            const std::vector <_CalcType_> &xcoeffs,
            const std::vector <_CalcType_> &ycoeffs
    ) {
        this->gain = gain;
        for (size_t i = 0; i < 2; ++i) {
            this->xcoeffs[i] = xcoeffs[i];
        }
        for (size_t i = 0; i < 2; ++i) {
            this->ycoeffs[i] = ycoeffs[i];
        }
    }
    void setFilterParameters(const FilterCoefficients &fc) {
        setFilterParameters(1.0 / fc.rGain, fc.xcoeffs, fc.ycoeffs);
    }

};
