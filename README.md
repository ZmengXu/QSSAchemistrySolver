# QSSAchemistrySolver
The quasi-steady-state approximation chemistry solver for OpenFOAM 4 and 7

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

