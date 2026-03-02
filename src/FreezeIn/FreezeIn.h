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

//Boost C++ library
#include <boost/math/special_functions/bessel.hpp>/*provides bessel-K function*/
#include <boost/math/quadrature/gauss.hpp>/*provides Gauss-Legendre quadrature*/
#include <boost/math/quadrature/exp_sinh.hpp>/*provides exp_sinh quadrature*/

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
long double M2_ffchichi(long double s, long double mchi, long double mf,long double gD, long double Nf, long double qH,long double ma) {

    //Axial (Af) piece of Aprime f f couplings
    long double Af = 0.5L*gD*qH;

    //Axial (Ac) piece of Aprime chi chi couplings
    long double Ac = 0.5L*gD*1.0L;

    return (16.0L * pow(Af, 2) * pow(Ac, 2) * pow(mf, 2) * pow(mchi, 2)) / pow(ma, 4);
}

/****************************************/
/* Collision terms for number densities */
/****************************************/

//Number-density collision term for f f -> Aprime/Z -> Chi Chi
long double CollisionNum_ffchichi(long double T, long double mchi,long double mf, long double gD,long double Nf, long double qH,long double ma, long double LambdaQCD) {

    if ( ( Nf == 1.0L ) || ( (Nf == 3.0L) && (T > LambdaQCD) ) ) {

        auto integrand_s = [=] (long double s) {
            return M2_ffchichi(s, mchi, mf, gD, Nf, qH, ma) * sqrt(1.0L - 4.0L*mchi*mchi/s) * sqrt(1.0L - 4.0L*mf*mf/s) * sqrt(s) * boost::math::cyl_bessel_k(1, sqrt(s)/T);
        };
        
        return (T/(pow(8.0L*M_PI, 2)*pow(2.0L*M_PI, 3))) *
               exp_sinh<long double>().integrate(integrand_s, max(4.0L*mf*mf, 4.0L*mchi*mchi), INFINITY);
    }
    else { return 0.0L; }
}

//Sum of all number-density collision terms for portal freeze-in
//qh1 for leptons, thetaD for quarks.
long double CollisionNum_chi(long double T, long double mchi, long double gD, long double qh1, long double thetaD, long double ma, long double anom_mass, long double LambdaQCD) {

    long double result = CollisionNum_ffchichi(T, mchi, Me, gD, 1.0L, qh1, ma, LambdaQCD) + /*e*/

                         CollisionNum_ffchichi(T, mchi, Mmu, gD, 1.0L, qh1, ma, LambdaQCD) + /*mu*/

                         CollisionNum_ffchichi(T, mchi, Mta, gD, 1.0L, qh1, ma, LambdaQCD); /*ta*/

    if (anom_mass != 0.0) {
        result = result  + CollisionNum_ffchichi(T, mchi, anom_mass, gD, 1.0L, qh1, ma, LambdaQCD) /*E*/;
    }

    return result;
}
// CollisionNum_ffchichi(T, mchi, Mu, gD, 3.0L, thetaD, ma, LambdaQCD) + /*u*/
// CollisionNum_ffchichi(T, mchi, Mc, gD, 3.0L, thetaD, ma, LambdaQCD) + /*c*/
// CollisionNum_ffchichi(T, mchi, Mt, gD, 3.0L, thetaD, ma, LambdaQCD) + /*t*/
// CollisionNum_ffchichi(T, mchi, Md, gD, 3.0L, thetaD, ma, LambdaQCD) + /*d*/
// CollisionNum_ffchichi(T, mchi, Ms, gD, 3.0L, thetaD, ma, LambdaQCD) + /*s*/
// CollisionNum_ffchichi(T, mchi, Mb, gD, 3.0L, thetaD, ma, LambdaQCD) /*b*/;
// CollisionNum_ffchichi(T, mchi, anom_mass, gD, 3.0L, thetaD, ma, LambdaQCD) /*U*/ 
// CollisionNum_ffchichi(T, mchi, anom_mass, gD, 3.0L, thetaD, ma, LambdaQCD) /*D*/

/*********************************************************/
/* Thermally-averaged cross-section for portal freeze-in */
/*********************************************************/

//Equilibrium number density for Chi
long double NumEq(long double T, long double m, int dof) {

    return (dof/(2.0L*M_PI*M_PI))*T*m*m*boost::math::cyl_bessel_k(2, m/T);

}

//Thermally-averaged cross section
long double SigmaV_chi(long double T, long double mchi, long double gD, long double qh1, long double thetaD, long double ma, long double anom_mass, long double LambdaQCD) {
    return CollisionNum_chi(T, mchi, gD, qh1, thetaD, ma, anom_mass,LambdaQCD) /
           pow(NumEq(T, mchi, 2), 2.0L);
}

/*****************************/
/* Freeze-in portal coupling */
/*****************************/

//Portal Yield for Chi
long double Yield_FreezeIn(long double mchi, long double gD, long double qh1, long double thetaD, long double ma, long double anom_mass, long double LambdaQCD, long double Trh) {

    auto integrand_T = [=] (long double T) {
        return HoverHbarVisible(T) *
               CollisionNum_chi(T, mchi, gD, qh1, thetaD, ma, anom_mass, LambdaQCD) /
               (gstarS(T)*sqrt(gstar(T))*pow(T, 6.0L));
    };
    return (135.0L*sqrt(10.0L)*MPl/(2.0L*pow(M_PI, 3.0L))) *
           gauss<long double, 701>().integrate(integrand_T, 0.0L, Trh);
}

//Portal coupling, gD, for freezing-in the required relic abundance
long double gD_FreezeIn(long double mchi, long double qh1, long double thetaD, long double ma, long double anom_mass, long double LambdaQCD, long double Trh) {
    if (Trh == 0.0L) {
        Trh = INFINITY;
    }
    return pow(
                4.37e-10L /
                (2.0L * mchi * Yield_FreezeIn(mchi, 1.0L, qh1, thetaD, ma, anom_mass, LambdaQCD, Trh)), 0.25L
               );
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
long double SigmaDDe(long double mchi, long double gD, long double qh1, long double ma) {

    long double Ae = 0.5L*gD*qh1;
    long double Ac = 0.5L*gD*1.0L;

    return (pow(Muchie(mchi), 2.0L) * pow(Ae, 2.0L) * pow(Ac, 2.0L)) / (M_PI * pow(ma, 4.0L))*pow(GeVinvtocm, 2.0L);
}

#endif
