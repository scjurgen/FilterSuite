//
// Created by Jürgen Schwietering on 10/31/15.
//
#include <gtest/gtest.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-generated-matchers.h>

using ::testing::ElementsAre;

#include "../FilterCoefficients.h"

const double epsilon = 1E-7;


TEST(FilterCoefficient, BESSEL_LowPass_1)
{
    FilterCoefficients f(PM_LP, FC_BESSEL, 1);
    f.computeFilter(0.01);
    EXPECT_NEAR(3.046874709125380054e-2,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}


TEST(FilterCoefficient, BESSEL_RecalcSameLowPass_1)
{
    FilterCoefficients f(PM_LP, FC_BESSEL, 1);
    f.computeFilter(0.01); // running twice should give same result
    f.computeFilter(0.01);
    EXPECT_NEAR(3.046874709125380054e-2,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}


TEST(FilterCoefficient, BESSEL_RecalcDifferentLowPass_1)
{
    FilterCoefficients f(PM_LP, FC_BESSEL, 1);
    f.computeFilter(0.1); // running again should give different result
    f.computeFilter(0.01);
    EXPECT_NEAR(3.046874709125380054e-2,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}


TEST(FilterCoefficient, BESSEL_HighPass_1)
{
    FilterCoefficients f(PM_HP, FC_BESSEL, 1);
    f.computeFilter(0.01);
    EXPECT_NEAR(9.695312529087461995e-1,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(-1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}

TEST(FilterCoefficient, BESSEL_BandPass_1)
{
    FilterCoefficients f(PM_BP, FC_BESSEL, 1);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(3.23764746428618934715e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1,epsilon),
            ::testing::DoubleNear(0,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39062505817492065852e-01,epsilon),
            ::testing::DoubleNear(1.93140991198042377697e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BESSEL_BandStop_1)
{
    FilterCoefficients f(PM_BS, FC_BESSEL, 1);
    f.computeFilter(0.01,0.02);

    EXPECT_NEAR(1.03142626604335507778e+00, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(9.99999999999999777955e-01,epsilon),
            ::testing::DoubleNear( -1.99210691371308601383e+00,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39062505817492065852e-01,epsilon),
            ::testing::DoubleNear(1.93140991198042377697e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, BUTTERWORTH_LowPass_1)
{
    FilterCoefficients f(PM_LP, FC_BUTTERWORTH, 1);
    f.computeFilter(0.01);
    EXPECT_NEAR(3.046874709125380054e-2,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}

TEST(FilterCoefficient, BUTTERWORTH_HighPass_1)
{
    FilterCoefficients f(PM_HP, FC_BUTTERWORTH, 1);
    f.computeFilter(0.01);
    EXPECT_NEAR(9.695312529087461995e-1,1.0/f.rGain,epsilon);
    EXPECT_NEAR(0.93906250581749239892,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(-1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}

TEST(FilterCoefficient, BUTTERWORTH_BandPass_1)
{
    FilterCoefficients f(PM_BP, FC_BUTTERWORTH, 1);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(3.23764746428618934715e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1,epsilon),
            ::testing::DoubleNear(0,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39062505817492065852e-01,epsilon),
            ::testing::DoubleNear(1.93140991198042377697e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, BUTTERWORTH_RecalcBandPass_1)
{
    FilterCoefficients f(PM_BP, FC_BUTTERWORTH, 1);
    f.computeFilter(0.1,0.2);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(3.23764746428618934715e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1,epsilon),
            ::testing::DoubleNear(0,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39062505817492065852e-01,epsilon),
            ::testing::DoubleNear(1.93140991198042377697e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BUTTERWORTH_BandStop_1)
{
    FilterCoefficients f(PM_BS, FC_BUTTERWORTH, 1);
    f.computeFilter(0.01,0.02);

    EXPECT_NEAR(1.03142626604335507778e+00, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(9.99999999999999777955e-01,epsilon),
            ::testing::DoubleNear( -1.99210691371308601383e+00,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39062505817492065852e-01,epsilon),
            ::testing::DoubleNear(1.93140991198042377697e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, CHEBYSHEV_LowPass_1)
{
    FilterCoefficients f(PM_LP, FC_CHEBYSHEV, 1);
    f.chebyshevRipple = -3;
    f.computeFilter(0.01);
    EXPECT_NEAR(3.27450486715411059890e+01,f.rGain,epsilon);
    EXPECT_NEAR(9.38922063605353240945e-01,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}

TEST(FilterCoefficient, CHEBYSHEV_HighPass_1)
{
    FilterCoefficients f(PM_HP, FC_CHEBYSHEV, 1);
    f.chebyshevRipple = -3;
    f.computeFilter(0.01);
    EXPECT_NEAR(1.03135173378584599213e+00,f.rGain,epsilon);
    EXPECT_NEAR(9.39202635223656989716e-01,f.ycoeffs[0],epsilon);
    EXPECT_NEAR(-1.0, f.xcoeffs[0],epsilon);
    EXPECT_NEAR(1.0, f.xcoeffs[1],epsilon);
}

TEST(FilterCoefficient, CHEBYSHEV_BandPass_1)
{
    FilterCoefficients f(PM_BP, FC_CHEBYSHEV, 1);
    f.chebyshevRipple = -3;
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(3.23040850703861082138e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1,epsilon),
            ::testing::DoubleNear(0,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.38922063605353129923e-01,epsilon),
            ::testing::DoubleNear(1.93127002402953418247e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, CHEBYSHEV_BandStop_1)
{
    FilterCoefficients f(PM_BS, FC_CHEBYSHEV, 1);
    f.chebyshevRipple = -3;
    f.computeFilter(0.01,0.02);

    EXPECT_NEAR(1.03135173378584443782e+00, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(9.99999999999999777955e-01,epsilon),
            ::testing::DoubleNear( -1.99210691371308601383e+00,epsilon),
            ::testing::DoubleNear(1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.39202635223656878694e-01,epsilon),
            ::testing::DoubleNear(1.93154948835984141553e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}



TEST(FilterCoefficient, RESONATOR_BandPass_1)
{
    FilterCoefficients f(PM_BP, FC_RESONATOR);
    f.qfactor = 10;
    f.computeFilter(0.01);
    EXPECT_NEAR(3.19310933380643348301e+02, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1,epsilon),
            ::testing::DoubleNear( 0,epsilon),
            ::testing::DoubleNear( 1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.93736512624778023373e-01,epsilon),
            ::testing::DoubleNear(1.98980232904289611184e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, RESONATOR_BandStop_1)
{
    FilterCoefficients f(PM_BS, FC_RESONATOR);
    f.qfactor = 10;
    f.computeFilter(0.01);
    EXPECT_NEAR(1.00314158231789507525e+00, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1,epsilon),
            ::testing::DoubleNear( -1.99605345685654311794e+00,epsilon),
            ::testing::DoubleNear( 1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.93736512624778023373e-01,epsilon),
            ::testing::DoubleNear(1.98980232904289611184e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, RESONATOR_AllPass_1)
{
    FilterCoefficients f(PM_AP, FC_RESONATOR);
    f.qfactor = 10;
    f.computeFilter(0.01);
    EXPECT_NEAR(1.00630296592270163103e+00, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.00630296592270518374e+00,epsilon),
            ::testing::DoubleNear( -2.00234398531577317826e+00,epsilon),
            ::testing::DoubleNear( 1,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.93736512624778023373e-01,epsilon),
            ::testing::DoubleNear(1.98980232904289611184e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BESSEL_LowPass_2)
{
    FilterCoefficients f(PM_LP, FC_BESSEL, 2);
    f.computeFilter(0.01);
    EXPECT_NEAR(6.70115907653785939146e+02, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 2.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-8.70683455111282866845e-01,epsilon),
            ::testing::DoubleNear(1.86471433849382361991e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, BUTTERWORTH_LowPass_2)
{
    FilterCoefficients f(PM_LP, FC_BUTTERWORTH, 2);
    f.computeFilter(0.01);
    EXPECT_NEAR(1.05854624078793835906e+03, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 2.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );
    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.14975834801433740573e-01,epsilon),
            ::testing::DoubleNear(1.91119706742607320393e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}



TEST(FilterCoefficient, BESSEL_BandPass_2)
{
    FilterCoefficients f(PM_BP, FC_BESSEL, 2);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(6.64373891939287318564e+02, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear(-2.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );

    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-8.70683455111283644001e-01,epsilon),
            ::testing::DoubleNear(3.59184969348952121138e+00,epsilon),
            ::testing::DoubleNear(-5.57068649457617137699e+00,epsilon),
            ::testing::DoubleNear(3.84946207661982242598e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}

TEST(FilterCoefficient, BUTTERWORTH_BandPass_2)
{
    FilterCoefficients f(PM_BP, FC_BUTTERWORTH, 2);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(1.05814276719598819909e+03, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear(-2.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );

    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-9.14975834801433629551e-01,epsilon),
            ::testing::DoubleNear(3.72638413213216646014e+00,epsilon),
            ::testing::DoubleNear(-5.70722925020648208516e+00,epsilon),
            ::testing::DoubleNear(3.89576135945691426343e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BUTTERWORTH_BandPass_2_B)
{
    FilterCoefficients f(PM_BP, FC_BUTTERWORTH, 2);
    f.computeFilter(0.1,0.2);
    EXPECT_NEAR(1.48233382093163132964e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear(-2.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );

    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-4.12801598096188548936e-01,epsilon),
            ::testing::DoubleNear(1.21665163551553057175e+00,epsilon),
            ::testing::DoubleNear(-2.11920239714428237932e+00,epsilon),
            ::testing::DoubleNear(1.94246877654788363543e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BUTTERWORTH_CopyBandPass_2_B)
{
    FilterCoefficients g(PM_BP, FC_BUTTERWORTH, 2);

    FilterCoefficients f(g);
    f.computeFilter(0.1,0.2);
    EXPECT_NEAR(1.48233382093163132964e+01, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear( 1.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear(-2.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );

    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-4.12801598096188548936e-01,epsilon),
            ::testing::DoubleNear(1.21665163551553057175e+00,epsilon),
            ::testing::DoubleNear(-2.11920239714428237932e+00,epsilon),
            ::testing::DoubleNear(1.94246877654788363543e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}


TEST(FilterCoefficient, BUTTERWORTH_3_B)
{
    FilterCoefficients g(PM_BP, FC_BUTTERWORTH, 3);

    FilterCoefficients f(g);
    f.computeFilter(0.01,0.02);
    EXPECT_NEAR(3.43090816488102718722e+04, f.rGain, epsilon);

    EXPECT_THAT(f.xcoeffs, ElementsAre(
            ::testing::DoubleNear(-1.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 3.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear(-3.0,epsilon),
            ::testing::DoubleNear( 0.0,epsilon),
            ::testing::DoubleNear( 1.0,epsilon)
    ) );

    EXPECT_THAT(f.ycoeffs, ElementsAre(
            ::testing::DoubleNear(-8.81893130592485308128e-01,epsilon),
            ::testing::DoubleNear(5.38084271928557456022e+00,epsilon),
            ::testing::DoubleNear(-1.37025732692692088222e+01,epsilon),
            ::testing::DoubleNear(1.86413712855399040791e+01,epsilon),
            ::testing::DoubleNear(-1.42889215565185185852e+01,epsilon),
            ::testing::DoubleNear(5.85117348976046969256e+00,epsilon),
            ::testing::DoubleNear(-1.0,epsilon)
    ) );
}
