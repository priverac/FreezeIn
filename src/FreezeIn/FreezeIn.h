#ifndef FREEZEIN_H
#define FREEZEIN_H

/********************************/
/* Standard and Boost Libraries */
/********************************/

//Standard libraries
#include <cmath>/*provides pow, sqrt, ...*/
#include <iostream>/*standard i/o library*/
#include <fstream>/*to read/write file*/
#include <cstdio>/*provides printf*/
#include <string>/*to use string data type*/
#include <vector>/*provides std::vector*/
#include <sstream>/*provides istringstream*/
#include <stdexcept>/*provides runtime_error*/
#include <map>/*provides std::map*/

//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/
#include <boost/math/quadrature/tanh_sinh.hpp>/*provides tanh_sinh quadrature*/

//Namespaces
using namespace std;
using namespace boost::math;
using namespace boost::math::quadrature;

/*********************************************************/
/* Masses, widths, coupling constants, and mixing angles */
/*********************************************************/

//Masses in GeV
#define Me 0.51099895e-3L
#define Mmu 105.6583755e-3L
#define Mta 1.77686L
#define Mu 2.16e-3L
#define Mc 1.27L
#define Mt 172.69L
#define Md 4.67e-3L
#define Ms 93.4e-3L
#define Mb 4.18L
#define Mpip 139.57039e-3L
#define MKp 493.677e-3L
#define MZ 91.1876L
#define MW 80.379L
#define MPl 2.435323077e18L

//Widths in GeV
#define WZ 2.4952L

//Alpha strong at mZ
#define alphaS 0.1179L

//Fine structure constant
#define alphaEM 7.2973525664e-3L

//Weinberg angle
#define sW2 0.23121L /*Sin-squared Weinberg angle*/
const long double sW = sqrt(sW2); /*Sin(ThetaW)*/
const long double cW = sqrt(1.0L - sW2); /*Cos(ThetaW)*/
const long double tW = sqrt(sW2/(1.0L - sW2)); /*Tan(ThetaW)*/
const long double s2W = 2.0L*sqrt(sW2 - sW2*sW2); /*Sin(2 ThetaW)*/

//Conversions
#define GeVinvtocm 1.97326937e-14 /*Inverse GeV in cm*/

/************************************************************************/
/* Useful fuctions: linear interpolation, finite-difference derivatives */
/************************************************************************/

//Linear interpolation function
long double interp(long double x, const vector<long double> &xData,
              const vector<long double> &yData, bool extrapolate) {

    if (xData.size() < 2 || yData.size() < 2) {
        throw runtime_error("interp: input arrays are too small. Did Read_gstar() load correctly?");
    }

    //Check if x is increasing/decreasing
    bool increasing = xData[1] > xData[0];

    //Perform binary search to find the interval for interpolation
    int mid, low = 0, high = xData.size() - 1;
    while (high - low > 1) {
        
        mid = (low + high)/2;

        if (
                (increasing && x >= xData[mid]) ||
                (!increasing && x <= xData[mid])
           ) {
            low = mid;
        }
        else {
            high = mid;
        }
    }

    //Points on either side (unless beyond ends)
    long double xL = xData[low], yL = yData[low];
    long double xR = xData[low+1], yR = yData[low+1];

    if ( !extrapolate ) {//if beyond ends of array and not extrapolating
       if ( (increasing && x < xL) || (!increasing && x > xL) ) yR = yL;
       if ( (increasing && x > xR) || (!increasing && x < xR) ) yL = yR;
    }

    long double dydx = (yR - yL)/(xR - xL);//gradient
    
    return yL + dydx * ( x - xL );//linear interpolation
}

//Array of slopes using Finite-difference formulas
vector<long double> slopearray(const vector<long double>& xData,
                          const vector<long double>& yData) {

    long double x, y, x1, y1, x2, y2;
    vector<long double> dydxData;

    for (size_t i = 0; i < xData.size(); i++) {

        x = xData[i]; y = yData[i];
        if ( i == 0 ) {
            x1 = xData[i+1]; y1 = yData[i+1]; x2 = xData[i+2]; y2 = yData[i+2];
        }
        else if ( i == xData.size() - 1 ) {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i-2]; y2 = yData[i-2];
        }
        else {
            x1 = xData[i-1]; y1 = yData[i-1]; x2 = xData[i+1]; y2 = yData[i+1];
        }

        dydxData.push_back(y*(2.0L*x - (x1 + x2))/(x - x1)/(x - x2) +
                           y1*(x - x2)/(x1 - x)/(x1 - x2) +
                           y2*(x - x1)/(x2 - x)/(x2 - x1));
    }

    return dydxData;
}

/***********************************************************/
/* g*(S): Effective number of degrees of freedom in the SM */
/***********************************************************/

//Define arrays for Temperature T and gstar(S) in the SM
vector<long double> Tvec;

vector<long double> gstarvec;
vector<long double> dlngstardlnTvec;

vector<long double> gstarSvec;
vector<long double> dlngstarSdlnTvec;

//Read gstar(S) data from various .tab files in the gstar folder, and compute
//their finite-difference derivatives as a function of T
void Read_gstar(const string& choice, const string& gstarpath) {
    
    long double r1, r2, r3;
    vector <long double> tempvec;

    //Clear global vectors before filling them
    Tvec.clear(), Tvec.shrink_to_fit();
    gstarvec.clear(), gstarvec.shrink_to_fit();
    gstarSvec.clear(), gstarSvec.shrink_to_fit();
    dlngstardlnTvec.clear(), dlngstardlnTvec.shrink_to_fit();
    dlngstarSdlnTvec.clear(), dlngstarSdlnTvec.shrink_to_fit();

    //Open a file corresponding to the input choice
    string filename;
    if (choice == "standard") { filename = "gstar/std.tab"; }
    else if (choice == "HP_A") { filename = "gstar/HP_A.tab"; }
    else if (choice == "HP_B") { filename = "gstar/HP_B.tab"; }
    else if (choice == "HP_B2") { filename = "gstar/HP_B2.tab"; }
    else if (choice == "HP_B3") { filename = "gstar/HP_B3.tab"; }
    else if (choice == "HP_C") { filename = "gstar/HP_C.tab"; }

    //Accept either:
    //1) gstarpath = project root (expects gstar/std.tab), or
    //2) gstarpath = direct gstar folder (expects std.tab)
    string filename_root = gstarpath + "/" + filename;
    string basename = filename.substr(filename.find_last_of('/') + 1);
    string filename_direct = gstarpath + "/" + basename;

    ifstream file(filename_root);
    string loaded_filename = filename_root;
    if (!file.is_open()) {
        file.clear();
        file.open(filename_direct);
        loaded_filename = filename_direct;
    }
    //If a file is open, read line-by-line to extract three values (r1, r2, r3)
    //from each line and store them in the corresponding global vectors
    if (file.is_open()) {
        string line;
        while (getline(file, line)) {
            if (line.empty() || line[0] == '#') continue;//Skip lines with '#'
            istringstream iss(line);
            iss >> r1 >> r2 >> r3;
            Tvec.push_back(r1);
            gstarvec.push_back(r3);
            gstarSvec.push_back(r2);
        }
        file.close();
    }
    else {
        throw runtime_error("Read_gstar: unable to open file: " + filename_root + " or " + filename_direct);
    }

    if (Tvec.size() < 3 || gstarvec.size() < 3 || gstarSvec.size() < 3) {
        throw runtime_error("Read_gstar: loaded fewer than 3 data rows from: " + loaded_filename);
    }

    tempvec = slopearray(Tvec, gstarvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstardlnTvec.push_back((Tvec[i]/gstarvec[i])*tempvec[i]);
    }
    tempvec = slopearray(Tvec, gstarSvec);
    for (size_t i = 0; i < Tvec.size(); i++) {
        dlngstarSdlnTvec.push_back((Tvec[i]/gstarSvec[i])*tempvec[i]);
    }

    return;
}

//g*
long double gstar(long double T) {
    return interp(T, Tvec, gstarvec, false);
}

//g*S
long double gstarS(long double T) {
    return interp(T, Tvec, gstarSvec, false);
}

//dlng*S/dlnT
long double dlngstarSdlnT(long double T) {
    return interp(T, Tvec, dlngstarSdlnTvec, false);
}

//dlng*/dlnT
long double dlngstardlnT(long double T) {
    return interp(T, Tvec, dlngstardlnTvec, false);
}

/***************************************************************************/
/* Energy density, comoving entropy, and Hubble rate in the Visible sector */
/***************************************************************************/

//Rho Visible
long double RhoVisible(long double T) {
    return (M_PI*M_PI/30.0L)*gstar(T)*pow(T, 4.0L);
}

//Comoving entropy of the Visible sector
long double EntropyVisible(long double T) {
    return (2.0L*M_PI*M_PI/45.0L)*gstarS(T)*T*T*T;
}

//Hubble rate
long double Hubble(long double T) {
    return sqrt(M_PI*M_PI*gstar(T)/90.0L)*T*T/MPl;
}

//(H / Hbar) to account for varying gstarS only in the Visible sector
long double HoverHbarVisible(long double T) {
    return (1.0L + (1.0L/3.0L)*dlngstarSdlnT(T));
}

/******************************************/
/* Fully averaged squared Matrix elements */
/******************************************/

//Fully averaged matrix element squared for f f -> Aprime -> chi chi
long double M2_ffchichi(long double s, long double mchi, long double mf,long double vD, long double Nf, long double qH,long double tb, long double thetaD) {
 
    //Axial (Af) piece of Aprime f f couplings
    long double Af = 0.5L*qH*thetaD;
 
    //Axial (Ac) piece of Aprime chi chi couplings
    long double Ac = 0.5L*1.0L;
 
    return (16.0L * pow(Af, 2) * pow(Ac, 2) * pow(mf, 2) * pow(mchi, 2)) / pow(vD, 4);
}

/****************************************/
/* Collision terms for number densities */
/****************************************/

//Number-density collision term for f f -> Aprime/Z -> Chi Chi
long double CollisionNum_ffchichi(long double T, long double mchi,long double mf, long double vD,long double Nf, long double qH,long double tb, long double LambdaQCD, long double thetaD) {
 
    if ( ( Nf == 1.0L ) || ( (Nf == 3.0L) && (T > LambdaQCD) ) ) {
 
        auto integrand_s = [=] (long double s) {
            return M2_ffchichi(s, mchi, mf, vD, Nf, qH, tb, thetaD) * sqrt(1.0L - 4.0L*mchi*mchi/s) * sqrt(1.0L - 4.0L*mf*mf/s) * sqrt(s) * boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };
        
        return (T/(pow(8.0L*M_PI, 2)*pow(2.0L*M_PI, 3))) *
               exp_sinh<long double>().integrate(integrand_s, max(4.0L*mf*mf, 4.0L*mchi*mchi), INFINITY);
    }
    else { return 0.0L; }
}

//Individual number-density collision terms for each fermion species.
//Returns a name->value map so each CollisionNum_ffchichi output is accessible.
map<string, long double> CollisionNum_chi_individual(long double T, long double mchi, long double vD, long double qh1, long double tb, long double anom_mass, long double LambdaQCD) {

    long double thetaL = (2.0L + pow(tb,2))/(1.0L + pow(tb,2)); /*leptons*/
    long double thetaQ = (1.0L)/(1.0L + pow(tb,2));             /*quarks*/

    map<string, long double> contribs;

    contribs["e"]  = CollisionNum_ffchichi(T, mchi, Me,  vD, 1.0L, qh1, tb, LambdaQCD, thetaL);
    contribs["mu"] = CollisionNum_ffchichi(T, mchi, Mmu, vD, 1.0L, qh1, tb, LambdaQCD, thetaL);
    contribs["ta"] = CollisionNum_ffchichi(T, mchi, Mta, vD, 1.0L, qh1, tb, LambdaQCD, thetaL);
    contribs["u"]  = CollisionNum_ffchichi(T, mchi, Mu,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);
    contribs["c"]  = CollisionNum_ffchichi(T, mchi, Mc,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);
    contribs["t"]  = CollisionNum_ffchichi(T, mchi, Mt,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);
    contribs["d"]  = CollisionNum_ffchichi(T, mchi, Md,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);
    contribs["s"]  = CollisionNum_ffchichi(T, mchi, Ms,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);
    contribs["b"]  = CollisionNum_ffchichi(T, mchi, Mb,  vD, 3.0L, qh1, tb, LambdaQCD, thetaQ);

    if (anom_mass != 0.0L) {
        contribs["E"] = CollisionNum_ffchichi(T, mchi, anom_mass, vD, 1.0L, qh1, tb, LambdaQCD, thetaL);
    }

    return contribs;
}

//Sum of all number-density collision terms for portal freeze-in
//qh1 for leptons.
long double CollisionNum_chi(long double T, long double mchi, long double vD, long double qh1, long double tb, long double anom_mass, long double LambdaQCD) {
 
    map<string, long double> contribs = CollisionNum_chi_individual(T, mchi, vD, qh1, tb, anom_mass, LambdaQCD);
 
    long double result = 0.0L;
    for (const auto& kv : contribs) result += kv.second;
 
    return result;
}

/*********************************************************/
/* Thermally-averaged cross-section for portal freeze-in */
/*********************************************************/

//Equilibrium number density for Chi
long double NumEq(long double T, long double m, int dof) {
 
    return (dof/(2.0L*M_PI*M_PI))*T*m*m*boost::math::cyl_bessel_k(2, m/T);
 
}
 
//Thermally-averaged cross section
long double SigmaV_chi(long double T, long double mchi, long double vD, long double qh1, long double tb, long double anom_mass, long double LambdaQCD) {
    return CollisionNum_chi(T, mchi, vD, qh1, tb, anom_mass, LambdaQCD) /
           pow(NumEq(T, mchi, 2), 2.0L);
}

/*****************************/
/* Freeze-in portal coupling */
/*****************************/

//Portal Yield for Chi
long double Yield_FreezeIn(long double mchi, long double vD, long double qh1, long double tb, long double anom_mass, long double LambdaQCD, long double Trh) {

    auto integrand_T = [=] (long double T) {
        return HoverHbarVisible(T) *
               CollisionNum_chi(T, mchi, vD, qh1, tb, anom_mass, LambdaQCD) /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };
    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L))) *
           gauss<long double, 701>().integrate(integrand_T, 0.0L, Trh);
}

//Portal coupling, gD, for freezing-in the required relic abundance
long double vD_FreezeIn(long double mchi, long double qh1, long double tb, long double anom_mass, long double LambdaQCD, long double Trh) {
    if (Trh == 0.0L) {
        Trh = INFINITY;
    }
    return pow(
                (2.0L * mchi * Yield_FreezeIn(mchi, 1.0L, qh1, tb, anom_mass, LambdaQCD, Trh)) /
                4.37e-10L, 0.25L
               );
}

//Running yield: integrate the freeze-in integrand from Tlow up to Thigh
long double Yield_FreezeIn_partial(long double mchi, long double vD, long double qh1, long double tb, long double anom_mass, long double LambdaQCD, long double Tlow, long double Thigh) {
 
    auto integrand_T = [=] (long double T) {
        return HoverHbarVisible(T) *
               CollisionNum_chi(T, mchi, vD, qh1, tb, anom_mass, LambdaQCD) /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };
    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L))) *
           gauss<long double, 701>().integrate(integrand_T, Tlow, Thigh);
}



/**********************************/
/* Direct detection cross section */
/**********************************/

//Reduced mass of Chi and electron
long double Muchie(long double mchi) {
    return mchi*Me/(mchi + Me);
}
 
//Direct detection cross section in squared-centimeter: \overline{\sigma}_e
//Isolated ma << 1 limit
long double SigmaDDe(long double mchi, long double vD, long double qh1, long double tb) {
 
    long double Ae = 0.5L*qh1*(2.0L+pow(tb, 2))/(1.0L+pow(tb, 2));
    long double Ac = 0.5L;
 
    return (2 * pow(Muchie(mchi), 2.0L) * pow(Ae, 2.0L) * pow(Ac, 2.0L)) / (M_PI * pow(vD, 4.0L))*pow(GeVinvtocm, 2.0L);
}

/*************************************************************************/
/* Vector-Vector (VV) case: f f -> Aprime -> chi chi                     */
/* Pure kinetic mixing, ALL charged SM fermions.                         */
/*                                                                       */
/* Matrix element is the Aprime piece of Eq. S4 in Bhattiprolu, McGehee  */
/* & Pierce (arXiv:2312.14152):                                          */
/*                                                                       */
/*   M2 = (32/3) Nf pi^2 alpha^2 kappa^2 Qf^2                            */
/*        * (1 + 2 mf^2/s) (1 + 2 mchi^2/s)                              */
/*                                                                       */
/* No Z exchange and no Aprime/Z interference, so neutrinos (Qf = 0) do  */
/* NOT contribute and are omitted. Also omitted, as before: pi+pi-,      */
/* K+K-, W+W-, and plasmon decays -- fermion channels only.              */
/*                                                                       */
/* The coupling is kappa = epsilon*sqrt(alpha'/alpha).                   */
/*                                                                       */
/* NOTE: unlike the AA/AV cases, Nf here is a genuine multiplier (the    */
/* color factor) as well as the LambdaQCD gate for quarks.               */
/*                                                                       */
/* Collision-term prefactor 16T/(4 pi)^5 matches the BMP convention for  */
/* this M2 -- do NOT swap in the AA/AV prefactor T/((8 pi)^2 (2 pi)^3).  */
/*                                                                       */
/* REPLACES the electron-only VV section (M2_eechichi_VV and the old     */
/* CollisionNum_chi_VV / yield / gD_FreezeIn functions).                 */
/*************************************************************************/

//Fully averaged matrix element squared for f f -> Aprime -> chi chi (VV)
long double M2_ffchichi_VV(long double s, long double mchi, long double mf,
                           long double kappa, long double Nf,
                           long double Qf) {

    return (32.0L/3.0L)*Nf*M_PI*M_PI*alphaEM*alphaEM*kappa*kappa*Qf*Qf*
           (1.0L + 2.0L*mf*mf/s)*(1.0L + 2.0L*mchi*mchi/s);
}

//Number-density collision term for f f -> Aprime -> chi chi (VV case)
long double CollisionNum_ffchichi_VV(long double T, long double mchi,
                                     long double mf, long double kappa,
                                     long double Nf, long double Qf,
                                     long double LambdaQCD) {

    if ( ( Nf == 1.0L ) || ( (Nf == 3.0L) && (T > LambdaQCD) ) ) {

        auto integrand_s = [=] (long double s) {
            return M2_ffchichi_VV(s, mchi, mf, kappa, Nf, Qf) *
                   sqrt(1.0L - 4.0L*mchi*mchi/s) *
                   sqrt(1.0L - 4.0L*mf*mf/s) *
                   sqrt(s) *
                   boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };

        return (16.0L*T/pow(4.0L*M_PI, 5.0L)) *
               exp_sinh<long double>().integrate(integrand_s,
                                            max(4.0L*mf*mf, 4.0L*mchi*mchi),
                                            INFINITY);
    }
    else { return 0.0L; }
}

//Individual number-density collision terms for each fermion species (VV).
//Returns a name->value map, mirroring CollisionNum_chi_individual.
map<string, long double> CollisionNum_chi_individual_VV(long double T,
        long double mchi, long double kappa, long double LambdaQCD) {

    map<string, long double> contribs;

    contribs["e"]  = CollisionNum_ffchichi_VV(T, mchi, Me,  kappa,
                                              1.0L, -1.0L, LambdaQCD);
    contribs["mu"] = CollisionNum_ffchichi_VV(T, mchi, Mmu, kappa,
                                              1.0L, -1.0L, LambdaQCD);
    contribs["ta"] = CollisionNum_ffchichi_VV(T, mchi, Mta, kappa,
                                              1.0L, -1.0L, LambdaQCD);
    contribs["u"]  = CollisionNum_ffchichi_VV(T, mchi, Mu,  kappa,
                                              3.0L, 2.0L/3.0L, LambdaQCD);
    contribs["c"]  = CollisionNum_ffchichi_VV(T, mchi, Mc,  kappa,
                                              3.0L, 2.0L/3.0L, LambdaQCD);
    contribs["t"]  = CollisionNum_ffchichi_VV(T, mchi, Mt,  kappa,
                                              3.0L, 2.0L/3.0L, LambdaQCD);
    contribs["d"]  = CollisionNum_ffchichi_VV(T, mchi, Md,  kappa,
                                              3.0L, -1.0L/3.0L, LambdaQCD);
    contribs["s"]  = CollisionNum_ffchichi_VV(T, mchi, Ms,  kappa,
                                              3.0L, -1.0L/3.0L, LambdaQCD);
    contribs["b"]  = CollisionNum_ffchichi_VV(T, mchi, Mb,  kappa,
                                              3.0L, -1.0L/3.0L, LambdaQCD);

    return contribs;
}

//Sum of all number-density collision terms for portal freeze-in (VV case)
long double CollisionNum_chi_VV(long double T, long double mchi,
                                long double kappa, long double LambdaQCD) {

    map<string, long double> contribs =
        CollisionNum_chi_individual_VV(T, mchi, kappa, LambdaQCD);

    long double result = 0.0L;
    for (const auto& kv : contribs) result += kv.second;

    return result;
}

//Running yield (VV case), integrated in u = ln(T).
//Guards: Tlow = 0 is floored at max(mchi, Me)/50. The electron is still
//the lightest contributing channel (neutrinos have Qf = 0 and drop out),
//so the lightest threshold is 2*max(mchi, Me) and the rate is
//Boltzmann-dead below this. Thigh = infinity is capped at 1e6 GeV, where
//the ~1/T^2 tail contributes at the ~1e-6 relative level.
long double Yield_FreezeIn_partial_VV(long double mchi, long double kappa,
                                      long double LambdaQCD,
                                      long double Tlow, long double Thigh) {

    long double Tmin = max(mchi, (long double)Me)/50.0L;
    if (Tlow < Tmin) Tlow = Tmin;

    if (isinf(Thigh)) Thigh = 1.0e6L;

    if (Thigh <= Tlow) return 0.0L;

    auto integrand_u = [=] (long double u) {
        long double T = exp(u);
        return T * HoverHbarVisible(T) *
               CollisionNum_chi_VV(T, mchi, kappa, LambdaQCD) /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };
    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L))) *
           gauss<long double, 701>().integrate(integrand_u,
                                               log(Tlow), log(Thigh));
}

//Portal Yield for Chi (VV case): full integral from T = 0 up to Trh
long double Yield_FreezeIn_VV(long double mchi, long double kappa,
                              long double LambdaQCD, long double Trh) {
    return Yield_FreezeIn_partial_VV(mchi, kappa, LambdaQCD, 0.0L, Trh);
}

//Portal coupling for freezing-in the required relic abundance (VV case).
//This is kappa = epsilon*sqrt(alpha'/alpha). Since M2 ~ kappa^2:
//    kappa = sqrt( 4.37e-10 / (2 mchi Y(kappa=1)) ).
//(The paper's kappa_FI = 1.94e-11 at mchi = 1 MeV additionally includes
//Z exchange, pi+pi-, K+K-, W+W-, and plasmon decays.)
long double gD_FreezeIn(long double mchi, long double LambdaQCD,
                        long double Trh) {
    if (Trh == 0.0L) {
        Trh = INFINITY;
    }
    return sqrt(
                4.37e-10L /
                (2.0L * mchi * Yield_FreezeIn_VV(mchi, 1.0L, LambdaQCD, Trh))
               );
}

/*************************************************************************/
/* Axial-Vector (AV) case: f f -> Aprime -> chi chi                      */
/*                                                                       */
/* M2 = 2 ( Achi^2 Vf^2 + 4 Af Achi Vf Vchi + Af^2 Vchi^2 )              */
/*                                                                       */
/* Vector couplings (species-dependent):                                 */
/*   V_lep  = (1/2) gD qH1 [ 1 + 2(-1/2 + 2 sW2)/(1 + tb^2) ]            */
/*   V_up   =       gD qH1 ( 1/2 - (4/3) sW2 )/(1 + tb^2)                */
/*   V_down =       gD qH1 (-1/2 + (2/3) sW2 )/(1 + tb^2)                */
/*   V_chi  = (1/2) gD                                                   */
/*                                                                       */
/* Axial couplings: carried over from the AA section with gD explicit,   */
/*   A_f   = (1/2) gD qH1 thetaD,   thetaL = (2+tb^2)/(1+tb^2)           */
/*                                  thetaQ =       1 /(1+tb^2)           */
/*   A_chi = (1/2) gD                                                    */
/*                                                                       */
/* Collision-term prefactor T/((8 pi)^2 (2 pi)^3) matches the AA-case    */
/* convention (this M2 is in the same normalization as M2_ffchichi) --   */
/* do NOT use the VV-case prefactor 16T/(4 pi)^5.                        */
/*                                                                       */
/* As in the AA case, Nf only gates the quark channels at LambdaQCD --   */
/* there is no N_c = 3 color multiplier.                                 */
/*                                                                       */
/* Nothing in the AA or VV sections is modified by this block.           */
/*************************************************************************/

//Fully averaged matrix element squared for f f -> Aprime -> chi chi (AV)
long double M2_ffchichi_AV(long double s, long double mchi, long double mf,
                           long double gD, long double Nf, long double qH,
                           long double tb, long double thetaD,
                           long double Vf) {
 
    //Axial pieces
    long double Af = 0.5L*gD*qH*thetaD;
    long double Ac = 0.5L*gD;
 
    //Vector piece of the Aprime chi chi coupling
    long double Vc = 0.5L*gD;
 
    //The two coupling combinations that appear
    long double AfVc2 = Af*Af*Vc*Vc;  /*Af^2 Vchi^2*/
    long double AcVf2 = Ac*Ac*Vf*Vf;  /*Achi^2 Vf^2*/
 
    long double mf2   = mf*mf;
    long double mchi2 = mchi*mchi;
 
    return (4.0L/3.0L)*(AfVc2 + AcVf2)
           + (8.0L/(3.0L*s))*( AfVc2*(mchi2 - 2.0L*mf2)
                             + AcVf2*(mf2 - 2.0L*mchi2) )
           - (32.0L*mf2*mchi2/(3.0L*s*s))*(AfVc2 + AcVf2);
}

//Number-density collision term for f f -> Aprime -> chi chi (AV case)
long double CollisionNum_ffchichi_AV(long double T, long double mchi,
                                     long double mf, long double gD,
                                     long double Nf, long double qH,
                                     long double tb, long double LambdaQCD,
                                     long double thetaD, long double Vf) {

    if ( ( Nf == 1.0L ) || ( (Nf == 3.0L) && (T > LambdaQCD) ) ) {

        auto integrand_s = [=] (long double s) {
            return M2_ffchichi_AV(s, mchi, mf, gD, Nf, qH, tb, thetaD, Vf) *
                   sqrt(1.0L - 4.0L*mchi*mchi/s) *
                   sqrt(1.0L - 4.0L*mf*mf/s) *
                   sqrt(s) *
                   boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };

        return (T/(pow(8.0L*M_PI, 2)*pow(2.0L*M_PI, 3))) *
               exp_sinh<long double>().integrate(integrand_s,
                                            max(4.0L*mf*mf, 4.0L*mchi*mchi),
                                            INFINITY);
    }
    else { return 0.0L; }
}

//Sum of all number-density collision terms for portal freeze-in (AV case)
long double CollisionNum_chi_AV(long double T, long double mchi,
                                long double gD, long double qh1,
                                long double tb, long double anom_mass,
                                long double LambdaQCD) {

    long double tb2 = pow(tb, 2);

    //Axial theta factors (as in the AA case)
    long double thetaL = (2.0L + tb2)/(1.0L + tb2); /*leptons*/
    long double thetaQ = (1.0L)/(1.0L + tb2);       /*quarks*/

    //Vector couplings
    long double Vlep  = 0.5L*gD*qh1*(1.0L + 2.0L*(-0.5L + 2.0L*sW2)/(1.0L + tb2));
    long double Vup   = gD*qh1*(0.5L - (4.0L/3.0L)*sW2)/(1.0L + tb2);
    long double Vdown = gD*qh1*(-0.5L + (2.0L/3.0L)*sW2)/(1.0L + tb2);

    long double result = 0.0L;

    result += CollisionNum_ffchichi_AV(T, mchi, Me,  gD, 1.0L, qh1, tb,
                                       LambdaQCD, thetaL, Vlep);  /*e*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mmu, gD, 1.0L, qh1, tb,
                                       LambdaQCD, thetaL, Vlep);  /*mu*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mta, gD, 1.0L, qh1, tb,
                                       LambdaQCD, thetaL, Vlep);  /*ta*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mu,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vup);   /*u*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mc,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vup);   /*c*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mt,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vup);   /*t*/
    result += CollisionNum_ffchichi_AV(T, mchi, Md,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vdown); /*d*/
    result += CollisionNum_ffchichi_AV(T, mchi, Ms,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vdown); /*s*/
    result += CollisionNum_ffchichi_AV(T, mchi, Mb,  gD, 3.0L, qh1, tb,
                                       LambdaQCD, thetaQ, Vdown); /*b*/

    if (anom_mass != 0.0L) {
        result += CollisionNum_ffchichi_AV(T, mchi, anom_mass, gD, 1.0L, qh1,
                                           tb, LambdaQCD, thetaL, Vlep); /*E*/
    }

    return result;
}

//Running yield (AV case), integrated in u = ln(T).
//Guards: Tlow = 0 is floored at max(mchi, Me)/50 (the lightest channel is
//the electron, opening at 2*max(mchi, Me), so the rate is Boltzmann-dead
//below this); Thigh = infinity is capped at 1e6 GeV, where the ~1/T^2 tail
//contributes at the ~1e-6 relative level.
long double Yield_FreezeIn_partial_AV(long double mchi, long double gD,
                                      long double qh1, long double tb,
                                      long double anom_mass,
                                      long double LambdaQCD,
                                      long double Tlow, long double Thigh) {

    long double Tmin = max(mchi, (long double)Me)/50.0L;
    if (Tlow < Tmin) Tlow = Tmin;

    if (isinf(Thigh)) Thigh = 1.0e6L;

    if (Thigh <= Tlow) return 0.0L;

    auto integrand_u = [=] (long double u) {
        long double T = exp(u);
        return T * HoverHbarVisible(T) *
               CollisionNum_chi_AV(T, mchi, gD, qh1, tb, anom_mass,
                                   LambdaQCD) /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };
    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L))) *
           gauss<long double, 701>().integrate(integrand_u,
                                               log(Tlow), log(Thigh));
}

//Portal Yield for Chi (AV case): full integral from T = 0 up to Trh
long double Yield_FreezeIn_AV(long double mchi, long double gD,
                              long double qh1, long double tb,
                              long double anom_mass, long double LambdaQCD,
                              long double Trh) {
    return Yield_FreezeIn_partial_AV(mchi, gD, qh1, tb, anom_mass, LambdaQCD,
                                     0.0L, Trh);
}

//Portal coupling, gD, for freezing-in the required relic abundance (AV).
//Every coupling in M2 carries one power of gD, so M2 ~ gD^4 and the yield
//scales as Y(gD) = gD^4 * Y(gD=1). The relic condition 2 mchi Y = 4.37e-10
//then gives
//    gD = ( 4.37e-10 / (2 mchi Y(gD=1)) )^{1/4}.
//Note the coupling sits in the NUMERATOR of the yield here, so this is the
//reciprocal of the vD_FreezeIn inversion (where Y ~ 1/vD^4).
long double gD_FreezeIn_AV(long double mchi, long double qh1, long double tb,
                           long double anom_mass, long double LambdaQCD,
                           long double Trh) {
    if (Trh == 0.0L) {
        Trh = INFINITY;
    }
    return pow(
                4.37e-10L /
                (2.0L * mchi * Yield_FreezeIn_AV(mchi, 1.0L, qh1, tb,
                                                 anom_mass, LambdaQCD, Trh)),
                0.25L
               );
}



#endif
