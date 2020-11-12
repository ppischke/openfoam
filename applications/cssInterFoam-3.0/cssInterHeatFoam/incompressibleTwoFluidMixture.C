/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "incompressibleTwoFluidMixture.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(incompressibleTwoFluidMixture, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//- Calculate and return the laminar viscosity
void Foam::incompressibleTwoFluidMixture::calcNu()
{
    nuModel1_->correct();
    nuModel2_->correct();

    const volScalarField limitedAlpha1
    (
        "limitedAlpha1",
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    // Average kinematic viscosity calculated from dynamic viscosity
    nu_ = mu()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
    
    // Inverse-average kinematic viscosity calculated from dynamic viscosity
    nuInv_ = muInv()/(limitedAlpha1*rho1_ + (scalar(1) - limitedAlpha1)*rho2_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::incompressibleTwoFluidMixture::incompressibleTwoFluidMixture
(
    const volVectorField& U,
    const surfaceScalarField& phi
)
:
    incompressibleTwoPhaseMixture
    (  
        U,
        phi
    ),
    
    nuInv_
    (
        IOobject
        (
            "nuInv",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedScalar("nuInv", dimensionSet(0, 2, -1, 0, 0), 0),
        calculatedFvPatchScalarField::typeName
    )
{
    incompressibleTwoFluidMixture::calcNu();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoFluidMixture::muInv() const
{
    const volScalarField limitedAlpha1
    (
        min(max(alpha1_, scalar(0)), scalar(1))
    );

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "muInv", 1.0/(limitedAlpha1/rho1_/nuModel1_->nu()
          + (scalar(1) - limitedAlpha1)/rho2_/nuModel2_->nu())
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoFluidMixture::muInvf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "muInvf", 1.0/(alpha1f/rho1_/fvc::interpolate(nuModel1_->nu())
          + (scalar(1) - alpha1f)/rho2_/fvc::interpolate(nuModel2_->nu()))
        )
    );
}


Foam::tmp<Foam::surfaceScalarField>
Foam::incompressibleTwoFluidMixture::nuInvf() const
{
    const surfaceScalarField alpha1f
    (
        min(max(fvc::interpolate(alpha1_), scalar(0)), scalar(1))
    );

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            "nuInvf", 1.0/(alpha1f/rho1_/fvc::interpolate(nuModel1_->nu())
            + (scalar(1) - alpha1f)/rho2_/fvc::interpolate(nuModel2_->nu()))
	        / (alpha1f*rho1_ + (scalar(1) - alpha1f)*rho2_)
        )
    );
}

// ec: i want the basic kinematic viscosities

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoFluidMixture::nu1() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "nu1",nuModel1_->nu()
        )
    );
}

Foam::tmp<Foam::volScalarField>
Foam::incompressibleTwoFluidMixture::nu2() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            "nu2",nuModel2_->nu()
        )
    );
}

// end ec


// ************************************************************************* //
