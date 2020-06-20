# QSSAchemistrySolver
The quasi-steady-state approximation chemistry solver for OpenFOAM 4 and 7
This repository is migrated from [adhiraj-dasgupta/unsupportedContribOF23x](https://github.com/adhiraj-dasgupta/unsupportedContribOF23x). For OpenFOAM-2.3, please refer to that repository.

This repository works for OpenFOAM-4.x (master branch) and OpenFOAM-7 (of7 branch).

Demo in case:
a) In $FOAM_CASE/system/controlDict, compile to get a QSS mechanism library,
   in shell typing the following command to compile mechanism

    ´gfortran -shared -fPIC -o libMechnism19sCH4.so ckwyp.f´

b) In $FOAM_CASE/system/controlDict
    libs
    (
        "./chemkin/libMechnism19sCH4.so"
        "libQSSchemistryModel.so"
    );

c) In $FOAM_CASE/constant/chemistryProperties
    chemistryType
    {
        chemistrySolver QSS;
        chemistryThermo psi;
    }
    QSSCoeffs
    {
        solver          seulex;
        absTol          1e-12;
        relTol          0.1;
    }

Here's the temperature profile results for 19 species QSS and 30 species mechanism.

![TemperatureTrace](OF_vs_CHEMKINII.eps)
