/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Code based on OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) Adhiraj Dasgupta & Shijie Xu
    \\/     M anipulation  |                     
-------------------------------------------------------------------------------
 License
     This file is a derivative work of OpenFOAM.
     OpenFOAM is free software: you can redistribute it and/or modify it
     under the terms of the GNU General Public License as published by
     the Free Software Foundation, either version 3 of the License, or
     (at your option) any later version.
     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
     for more details.
     You should have received a copy of the GNU General Public License
     along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::QSS

Description
    A solver for chemistry when the mechanism is available as a Fortran ckwyp
    routine. Reduced mechanisms like these contain QSS species which are not
    transported, but only occur in the chemical source terms. Since fewer
    species are transported, this can lead to a saving in computational
    effort.

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

SourceFiles
    QSS.C

\*---------------------------------------------------------------------------*/

#ifndef QSS_H
#define QSS_H

#include "chemistrySolver.H"
#include "ODESolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

        //- the external routine to compute the reaction rates
        extern "C"
        {
            void ckwyp_
            (
                double* p,
                double* T,
                double  [],
                int* ickwrk,
                double* rckwrk,
                double [] 
            );
        }

/*---------------------------------------------------------------------------*\
                       Class QSS Declaration
\*---------------------------------------------------------------------------*/

template<class ChemistryModel>
class QSS
:
    public chemistrySolver<ChemistryModel>
{
    // Private data

        dictionary coeffsDict_;
        autoPtr<ODESolver> odeSolver_;

        // Solver data
        mutable scalarField cTp_;
        mutable double* Ytmp0;
        mutable double* dYdt;


public:

    //- Runtime type information
    TypeName("QSS");


    // Constructors

        //- Construct from mesh and phase name
        QSS(const fvMesh& mesh, const word& phaseName);


    //- Destructor
    virtual ~QSS();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual void solve
        (
            scalarField& c,
            scalar& T,
            scalar& p,
            scalar& deltaT,
            scalar& subDeltaT
        ) const;

        void jacobian
        (
            const scalar t,
            const scalarField& c,
            scalarField& dcdt,
            scalarSquareMatrix& dfdc
        ) const;
        
        tmp<scalarField> omega
        (
            const scalarField& c,
            const scalar T,
            const scalar p
        ) const;
        
        template<class ThermoType>
        scalar omega
        (
            const Reaction<ThermoType>& r,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;

        virtual scalar omegaI
        (
            label iReaction,
            const scalarField& c,
            const scalar T,
            const scalar p,
            scalar& pf,
            scalar& cf,
            label& lRef,
            scalar& pr,
            scalar& cr,
            label& rRef
        ) const;
        
        virtual tmp<DimensionedField<scalar, volMesh> > calculateRR
        (
            const label reactionI,
            const label specieI    
        ) const;
        
        virtual tmp<volScalarField> tc() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "QSS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
