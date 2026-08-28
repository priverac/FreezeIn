//This C++ program exposes the functions in the FreezeIn.h library to python
//using a light-weight header-only pybind11 library (included in this
//repository)

/********************/
/* FreezeIn Library */
/********************/

#include "FreezeIn.h"

/******************/
/* Pybind Library */
/******************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

#ifndef GSTARPATH
#define GSTARPATH "/Users/priverac/Documents/uiuc/FreezeIn/gstar"
#endif

/************************************/
/* Exposing C++ functions to python */
/************************************/

PYBIND11_MODULE(FreezeIn, mod)
{
    //Read_gstar(choice, gstarpath)
    mod.def("Read_gstar", &Read_gstar, R"pbdoc(
    Read tabulated data for effective number of degrees of freedom from various
    .tab files in the gstar folder with three columns:
    {Temperature in GeV, g*S, g*}.

    Inputs
    ------

    choices: "standard": Gondolo-Gelmini (LambdaQCD = 150 MeV) (default)
             "HP_A": Hindmarsh-Philipsen equation of state A
             "HP_B": Hindmarsh-Philipsen equation of state B
             "HP_B2": Hindmarsh-Philipsen equation of state B2
             "HP_B3": Hindmarsh-Philipsen equation of state B3
             "HP_C": Hindmarsh-Philipsen equation of state C
             (These are taken directly from MicrOMEGAs package)
    
    gstarpath: By default, set to the path to the gstar folder provided with
               this package
    )pbdoc", py::arg("choice")="standard", py::arg("gstarpath")=GSTARPATH);
    
    //gstar(T)
    mod.def("gstar", &gstar, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The effective number of degrees of freedom for energy density g*

    (By default uses the standard Gondolo-Gelmini g*(T). To use other choices
    for g*: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //gstarS(T)
    mod.def("gstarS", &gstarS, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The effective number of degrees of freedom for entropy density g*S

    (By default uses the standard Gondolo-Gelmini g*S(T). To use other choices
    for g*S: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //dlngstarSdlnT(T)
    mod.def("dlngstarSdlnT", &dlngstarSdlnT, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    The derivative of log(g*S) with respect to log(T), where g*S is the
    effective number of degrees of freedom for entropy density

    (By default uses the standard Gondolo-Gelmini g*S(T). To use other choices
    for g*S: evaluate Read_gstar(choice); see documentation for the function
    Read_gstar for more details.)
    )pbdoc", py::arg("T"));
    
    //RhoVisible(T)
    mod.def("RhoVisible", &RhoVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Energy density in the visible sector
    )pbdoc", py::arg("T"));
    
    //EntropyVisible(T)
    mod.def("EntropyVisible", &EntropyVisible, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Comoving entropy in the visible sector
    )pbdoc", py::arg("T"));

    //Hubble(T)
    mod.def("Hubble", &Hubble, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV

    Returns
    -------

    Hubble rate
    )pbdoc", py::arg("T"));

    mod.def("CollisionNum_chi_individual", &CollisionNum_chi_individual,
        py::arg("T"), py::arg("mchi"), py::arg("vD"), py::arg("qh1"), py::arg("tb"), py::arg("anom_mass"), py::arg("LambdaQCD"));


    //SigmaV_chi(T, mchi, gD, qh1, tb, ma, anom_mass, LambdaQCD, thetaD)
    mod.def("SigmaV_chi", &SigmaV_chi, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV
    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    ma: dark photon mass
    anom_mass: anomalon mass scale. Set to 0 = no anomalons by default
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default
    thetaD: coupling parameter for theta_D

    Returns
    -------

    Thermally-averaged cross section for SM SMbar -> chi chibar process
    )pbdoc", py::arg("T"), py::arg("mchi"), py::arg("vD"), py::arg("qh1"), py::arg("tb"), py::arg("anom_mass"), py::arg("LambdaQCD"));

    //vD_FreezeIn(mchi, qh1, tb, ma, anom_mass, LambdaQCD, Trh, thetaD)
    mod.def("vD_FreezeIn", &vD_FreezeIn, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    ma: mass of dark photon in GeV
    anom_mass: anomalon mass scale. Set to 0 = no anomalons by default
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default
    Trh: instantaneous reheating temperature. Setting to 0.0 will set it to infinity
    thetaD: coupling parameter for theta_D

    Returns
    -------

    Portal coupling gD that reproduces the observed dark matter relic
    abundance for dark matter frozen-in via a light dark photon mediator
    )pbdoc", py::arg("mchi"), py::arg("qh1"), py::arg("tb"), py::arg("anom_mass"), py::arg("LambdaQCD"), py::arg("Trh"));

    //Yield_FreezeIn(mchi, gD, qh1, tb, ma, anom_mass, LambdaQCD, Trh, thetaD)
    mod.def("Yield_FreezeIn", &Yield_FreezeIn, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    ma: mass of dark photon in GeV
    anom_mass: anomalon mass scale.
    LambdaQCD: QCD confinement scale in GeV. Set to 0.15 GeV by default
    Trh: instantaneous reheating temperature. Setting to 0.0 will set it to infinity
    thetaD: coupling parameter for theta_D

    Returns
    -------

    Portal yield for freeze-in dark matter as a function of temperature

    )pbdoc", py::arg("mchi"), py::arg("vD"), py::arg("qh1"), py::arg("tb"), py::arg("anom_mass"), py::arg("LambdaQCD"), py::arg("Trh"));

    mod.def("Yield_FreezeIn_partial", &Yield_FreezeIn_partial,
        py::arg("mchi"), py::arg("vD"), py::arg("qh1"), py::arg("tb"),
        py::arg("anom_mass"), py::arg("LambdaQCD"),
        py::arg("Tlow"), py::arg("Thigh"));

    //SigmaDDe(mchi, gD, qh1, tb, ma)
    mod.def("SigmaDDe", &SigmaDDe, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    ma: mass of dark photon

    Returns
    -------

    Direct detection cross section (in cm^2) through the light dark photon
    mediator
    )pbdoc", py::arg("mchi"), py::arg("vD"), py::arg("qh1"), py::arg("tb"));

    
    /**********************************************************/
    /* Vector-Vector (VV) case: e+ e- -> Aprime -> chi chi    */
    /**********************************************************/
 
    //CollisionNum_chi_VV(T, mchi, kappa, LambdaQCD)
    mod.def("CollisionNum_chi_VV", &CollisionNum_chi_VV,
        py::arg("T"), py::arg("mchi"), py::arg("kappa"),
        py::arg("LambdaQCD"));

    //CollisionNum_chi_individual_VV(T, mchi, kappa, LambdaQCD)
    mod.def("CollisionNum_chi_individual_VV", &CollisionNum_chi_individual_VV,
        py::arg("T"), py::arg("mchi"), py::arg("kappa"),
        py::arg("LambdaQCD"));

    //Yield_FreezeIn_VV(mchi, kappa, LambdaQCD, Trh)
    mod.def("Yield_FreezeIn_VV", &Yield_FreezeIn_VV,
        py::arg("mchi"), py::arg("kappa"), py::arg("LambdaQCD"),
        py::arg("Trh"));

    //Yield_FreezeIn_partial_VV(mchi, kappa, LambdaQCD, Tlow, Thigh)
    mod.def("Yield_FreezeIn_partial_VV", &Yield_FreezeIn_partial_VV,
        py::arg("mchi"), py::arg("kappa"), py::arg("LambdaQCD"),
        py::arg("Tlow"), py::arg("Thigh"));

    //gD_FreezeIn(mchi, LambdaQCD, Trh)
    mod.def("gD_FreezeIn", &gD_FreezeIn,
        py::arg("mchi"), py::arg("LambdaQCD"), py::arg("Trh"));

    /**********************************************************/
    /* Axial-Vector (AV) case: f f -> Aprime -> chi chi       */
    /**********************************************************/

    //CollisionNum_chi_AV(T, mchi, gD, qh1, tb, anom_mass, LambdaQCD)
    mod.def("CollisionNum_chi_AV", &CollisionNum_chi_AV, R"pbdoc(
    Inputs
    ------

    T: Temperature in the visible sector in GeV
    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    anom_mass: anomalon mass scale. Set to 0 = no anomalons
    LambdaQCD: QCD confinement scale in GeV

    Returns
    -------

    Sum of the number-density collision terms for SM SMbar -> chi chibar
    in the Axial-Vector case (all fermion channels)
    )pbdoc", py::arg("T"), py::arg("mchi"), py::arg("gD"), py::arg("qh1"),
             py::arg("tb"), py::arg("anom_mass"), py::arg("LambdaQCD"));

    //Yield_FreezeIn_AV(mchi, gD, qh1, tb, anom_mass, LambdaQCD, Trh)
    mod.def("Yield_FreezeIn_AV", &Yield_FreezeIn_AV, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    anom_mass: anomalon mass scale. Set to 0 = no anomalons
    LambdaQCD: QCD confinement scale in GeV
    Trh: instantaneous reheating temperature. Setting to 0.0 will NOT be
         remapped here; pass the desired upper limit explicitly

    Returns
    -------

    Portal yield for freeze-in dark matter in the Axial-Vector case,
    integrated from T = 0 up to Trh
    )pbdoc", py::arg("mchi"), py::arg("gD"), py::arg("qh1"), py::arg("tb"),
             py::arg("anom_mass"), py::arg("LambdaQCD"), py::arg("Trh"));

    //Yield_FreezeIn_partial_AV(mchi, gD, qh1, tb, anom_mass, LambdaQCD,
    //                          Tlow, Thigh)
    mod.def("Yield_FreezeIn_partial_AV", &Yield_FreezeIn_partial_AV, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    gD: portal coupling
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    anom_mass: anomalon mass scale. Set to 0 = no anomalons
    LambdaQCD: QCD confinement scale in GeV
    Tlow: lower limit of the temperature integral in GeV
    Thigh: upper limit of the temperature integral in GeV

    Returns
    -------

    Running (cumulative) portal yield in the Axial-Vector case,
    integrated from Tlow up to Thigh
    )pbdoc", py::arg("mchi"), py::arg("gD"), py::arg("qh1"), py::arg("tb"),
             py::arg("anom_mass"), py::arg("LambdaQCD"),
             py::arg("Tlow"), py::arg("Thigh"));

    //gD_FreezeIn_AV(mchi, qh1, tb, anom_mass, LambdaQCD, Trh)
    mod.def("gD_FreezeIn_AV", &gD_FreezeIn_AV, R"pbdoc(
    Inputs
    ------

    mchi: mass of the dark matter in GeV
    qh1: lepton Higgs charge
    tb: ratio of Higgs vacuum expectation values (tan(beta))
    anom_mass: anomalon mass scale. Set to 0 = no anomalons
    LambdaQCD: QCD confinement scale in GeV
    Trh: instantaneous reheating temperature. Setting to 0.0 will set it
         to infinity

    Returns
    -------

    Portal coupling gD that reproduces the observed dark matter relic
    abundance in the Axial-Vector case
    )pbdoc", py::arg("mchi"), py::arg("qh1"), py::arg("tb"),
             py::arg("anom_mass"), py::arg("LambdaQCD"), py::arg("Trh"));

};