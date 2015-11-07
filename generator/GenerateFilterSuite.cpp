
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "../FilterCoefficients.h"



std::string ReplaceString(std::string subject
        , const std::string& search
        , const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
        subject.replace(pos, search.length(), replace);
        pos += replace.length();
    }
    return subject;
}


std::string createFromTemplate(const std::string templateFilename, const std::map<std::string, std::string> &substitutions)
{
    std::string templateString;
    std::ifstream in;
    in.open(templateFilename,std::ios_base::binary|std::ios_base::in);
    if (!in.is_open())
    {
        std::cerr << "could not open FilterStringTemplate.h, set working directory to []/generator" << std::endl;
        exit(ENOENT);
    }
    std::stringstream buffer;
    buffer << in.rdbuf();
    templateString = buffer.str();
    for (auto item:substitutions)
    {
        templateString = ReplaceString(templateString, item.first, item.second);
    }
    return templateString;
}




void GenerateFilter(bool testonly, std::ofstream &outfileInclude, std::ofstream &outfileTest, PASSMODE pm, FILTERCHARACTER filtercharacter, int order)
{

    FilterCoefficients filterCoefficients(pm,filtercharacter,order);
    filterCoefficients.qfactor = 10;
    filterCoefficients.chebyshevRipple = -3;
    filterCoefficients.computeFilter(0.01,0.02);
    // same filter different alpha, checks which coeffcients are invariant
    FilterCoefficients filterCoefficientsInvariant(pm,filtercharacter,order);
    filterCoefficientsInvariant.qfactor = 10;
    filterCoefficientsInvariant.chebyshevRipple = -3;
    filterCoefficientsInvariant.computeFilter(0.12,0.23);

    std::stringstream xchain;
    std::stringstream ychain;
    for (size_t i=0; i < filterCoefficients.ycoeffs.size()-1; ++i)
    {
        ychain << " + ";
        if (filterCoefficients.ycoeffs[i]!=filterCoefficientsInvariant.ycoeffs[i])
            ychain << "v[" << i << "] * " << "ycoeffs[" << i << "]";
        else
            ychain << "v[" << i << "] * " << filterCoefficients.ycoeffs[i];
    }
    for (size_t i=0; i < filterCoefficients.xcoeffs.size(); ++i)
    {
        if (filterCoefficients.xcoeffs[i] !=0) {
            if (i>0)
                xchain << " + ";
            if (filterCoefficients.xcoeffs[i]!=filterCoefficientsInvariant.xcoeffs[i])
                xchain << "v[" << i << "] * " << "xcoeffs[" << i << "]";
            else {
                if (filterCoefficients.xcoeffs[i] == -1)
                    xchain << "-v[" << i << "]";
                else
                if (filterCoefficients.xcoeffs[i] == 1)
                    xchain << "v[" << i << "]";
                else
                    xchain << filterCoefficients.xcoeffs[i] << "*v[" << i << "]";
            }
        }
    }

    std::map<std::string, std::string> substitutions;
    switch(pm)
    {
        case PM_LP:substitutions.insert(std::pair<std::string,std::string>("__PASSMODE__", "PM_LP")); break;
        case PM_HP:substitutions.insert(std::pair<std::string,std::string>("__PASSMODE__", "PM_HP")); break;
        case PM_BP:substitutions.insert(std::pair<std::string,std::string>("__PASSMODE__", "PM_BP")); break;
        case PM_BS:substitutions.insert(std::pair<std::string,std::string>("__PASSMODE__", "PM_BS")); break;
        case PM_AP:substitutions.insert(std::pair<std::string,std::string>("__PASSMODE__", "PM_AP")); break;
    }

    if (pm == PM_BP || pm==PM_BS)
    {
        substitutions.insert(std::pair<std::string,std::string>(
                "__SETALPHA__",
                "void setAlphas(double aLow, double aHigh){filterCoefficients.alpha1=aLow; filterCoefficiencts.alpha2=aHigh;}"));
    }
    else
    {
        substitutions.insert(std::pair<std::string,std::string>(
                "__SETALPHA__",
                "void setAlphas(double alpha){filterCoefficients.alphaLow=aLow; filterCoefficiencts.alphaHigh=aLow;}"));
    }

    switch(filtercharacter)
    {
        case FC_BUTTERWORTH:substitutions.insert(std::pair<std::string,std::string>("__CHARACTER__", "FC_BUTTERWORTH")); break;
        case FC_BESSEL:substitutions.insert(std::pair<std::string,std::string>("__CHARACTER__", "FC_BESSEL")); break;
        case FC_CHEBYSHEV:substitutions.insert(std::pair<std::string,std::string>("__CHARACTER__", "FC_CHEBYSHEV")); break;
        case FC_RESONATOR:substitutions.insert(std::pair<std::string,std::string>("__CHARACTER__", "FC_RESONATOR")); break;
        case FC_PROPORTIONAL_INTEGRAL:substitutions.insert(std::pair<std::string,std::string>("__CHARACTER__", "FC_PROPORTIONAL_INTEGRAL")); break;
    }


    substitutions.insert(std::pair<std::string,std::string>("__EXTRAPARAMETER__", ""));
    substitutions.insert(std::pair<std::string,std::string>("__EXTRAPARAMETERPASS__", ""));

    substitutions.insert(std::pair<std::string,std::string>("__ORDER__", std::to_string(order)));
    substitutions.insert(std::pair<std::string,std::string>("__RIPLE__", std::to_string(filterCoefficients.chebyshevRipple)));
    substitutions.insert(std::pair<std::string,std::string>("__QFACTOR__", std::to_string(filterCoefficients.qfactor)));
    substitutions.insert(std::pair<std::string,std::string>("__FILTERNAME__", filterCoefficients.getCurrentPureFilterName()));
    substitutions.insert(std::pair<std::string,std::string>("__DELAYELEMENTS__", std::to_string(filterCoefficients.xcoeffs.size())));
    substitutions.insert(std::pair<std::string,std::string>("__DELAYELEMENTSM1__", std::to_string(filterCoefficients.xcoeffs.size()-1)));
    substitutions.insert(std::pair<std::string,std::string>("__XCOUNT__", std::to_string(filterCoefficients.ycoeffs.size()-1)));
    substitutions.insert(std::pair<std::string,std::string>("__XCHAIN__", xchain.str()));
    substitutions.insert(std::pair<std::string,std::string>("__YCHAIN__", ychain.str()));
    if (!testonly)
        outfileInclude << createFromTemplate("FilterStringTemplate.h", substitutions);
    outfileTest << createFromTemplate("FilterTestStringTemplate.h", substitutions);
    if (!testonly) {
        std::cout << filterCoefficients.getCurrentPureFilterName() << std::endl;
    }
}


int main(void)
{
    std::ofstream outfileInclude;
    outfileInclude.open("../AllFilters.h", std::ios_base::out | std::ios_base::binary);
    outfileInclude << "/** This file is autogenerated, do not change manually\n*/\n";

    std::ofstream outfileTest;
    outfileTest.open("../tests/AllFilters_unittest.cpp", std::ios_base::out | std::ios_base::binary);
    outfileTest << "/** This file is autogenerated, do not change manually\n*/\n";
    outfileTest << "\n"
                       "#include <gtest/gtest.h>\n"
                       "#include <gmock/gmock-matchers.h>\n"
                       "#include <gmock/gmock-generated-matchers.h>\n"
                       "using ::testing::ElementsAre;\n"
                       "#include \"../FilterCoefficients.h\"\n"
                       "#include \"../FilterSuite.h\"\n"
                       "#include \"FftBinCompare.h\"\n"
                       "#include \"CalcBins.h\"\n\n";
    std::vector<PASSMODE> passmodes = {
        PM_LP,
        PM_HP,
        PM_BP,
        PM_BS};

    for (int order=1; order < 4; order++)
    {
        for (auto pm:passmodes)
        {
            GenerateFilter(false, outfileInclude, outfileTest, pm, FC_BESSEL, order);
            GenerateFilter(true, outfileInclude, outfileTest, pm, FC_BUTTERWORTH, order);
            GenerateFilter(true, outfileInclude, outfileTest, pm, FC_CHEBYSHEV, order);
        }
    }
    std::vector<PASSMODE> passmodesResonator = {
            PM_AP,
            PM_BP,
            PM_BS};
    for (auto pm:passmodesResonator)
    {
        GenerateFilter(false,outfileInclude, outfileTest, pm, FC_RESONATOR, 1);
    }

    return 0;
}
