/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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
    Foam::incompressible::laminarTwoFluid

Description
    Turbulence model for laminar two fluid incompressible flow.

SourceFiles
    laminarTwoFluid.C

\*---------------------------------------------------------------------------*/

#ifndef laminarTwoFluid_H
#define laminarTwoFluid_H

#include "turbulenceModel.H" //deleted
//#include "turbulentTransportModel.H" //added
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                           Class laminarTwoFluid Declaration
\*---------------------------------------------------------------------------*/

class laminarTwoFluid
:
    public turbulenceModel    
{

    // Private data

        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

    
protected:

    // Protected data

        //- Interfacial viscosity tensor
		volSymmTensorField nuTensor_;

        
public:

    //- Runtime type information
    TypeName("laminarTwoFluid");

    // Constructors

        //- Construct from components
        laminarTwoFluid
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<laminarTwoFluid> New
        (
            const volVectorField& U,
            const surfaceScalarField& phi,
            transportModel& transport,
            const word& turbulenceModelName = turbulenceModel::typeName
        );


    //- Destructor
    virtual ~laminarTwoFluid()
    {}


    // Member Functions

        //- Const access to the coefficients dictionary
        virtual const dictionary& coeffDict() const;

        //- Return the turbulence viscosity, i.e. 0 for laminar flow
        virtual tmp<volScalarField> nut() const;

        //- Return the effective viscosity, i.e. the laminar viscosity
        virtual tmp<volScalarField> nuEff() const;

        //- Return the effective viscosity, i.e. the laminar viscosity
        virtual const volSymmTensorField& nuTensor() const;
                
        //- Return the turbulence kinetic energy, i.e. 0 for laminar flow
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate,
        //  i.e. 0 for laminarTwoFluid flow
        virtual tmp<volScalarField> epsilon() const;

        //- Return the Reynolds stress tensor, i.e. 0 for laminar flow
        virtual tmp<volSymmTensorField> R() const;

        //- Return the effective stress tensor, i.e. the laminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Return the source term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevRhoReff
        (
            const volScalarField& rho,
            volVectorField& U
        ) const;

        //- Correct the laminarTwoFluid viscosity
        virtual void correct();

        //- Read turbulenceProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
