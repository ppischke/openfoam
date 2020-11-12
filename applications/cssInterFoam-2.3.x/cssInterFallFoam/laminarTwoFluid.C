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

\*---------------------------------------------------------------------------*/

#include "laminarTwoFluid.H"
#include "incompressibleTwoFluidMixture.H"
#include "Time.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(laminarTwoFluid, 0);
addToRunTimeSelectionTable(turbulenceModel, laminarTwoFluid, turbulenceModel);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

laminarTwoFluid::laminarTwoFluid
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName
)
:
    turbulenceModel(U, phi, transport, turbulenceModelName),
    
    deltaN_
    (
        "deltaN",
        SMALL/pow(average(U_.mesh().V()), 1.0/3.0)
    ),
        
    nuTensor_
    (
        IOobject
        (
            "nuTensor",
            U_.time().timeName(),
            U_.db()
        ),
        U_.mesh(),
        dimensionedSymmTensor("nuTensor", dimensionSet(0, 2, -1, 0, 0), symmTensor::zero),
        calculatedFvPatchScalarField::typeName
    )
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<laminarTwoFluid> laminarTwoFluid::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName
)
{
    return autoPtr<laminarTwoFluid>
    (
        new laminarTwoFluid(U, phi, transport, turbulenceModelName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const dictionary& laminarTwoFluid::coeffDict() const
{
    return dictionary::null;
}


tmp<volScalarField> laminarTwoFluid::nut() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "nut",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("nut", nu()().dimensions(), 0.0)
        )
    );
}


tmp<volScalarField> laminarTwoFluid::nuEff() const
{
    return tmp<volScalarField>(new volScalarField("nuEff", nu()));
}


const volSymmTensorField& laminarTwoFluid::nuTensor() const
{
    return nuTensor_;
}


tmp<volScalarField> laminarTwoFluid::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("k", sqr(U_.dimensions()), 0.0)
        )
    );
}


tmp<volScalarField> laminarTwoFluid::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "epsilon", sqr(U_.dimensions())/dimTime, 0.0
            )
        )
    );
}


tmp<volSymmTensorField> laminarTwoFluid::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor
            (
                "R", sqr(U_.dimensions()), symmTensor::zero
            )
        )
    );
}


tmp<volSymmTensorField> laminarTwoFluid::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nu()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> laminarTwoFluid::divDevReff(volVectorField& U) const
{
    return
    (    
      - fvm::laplacian(nuTensor(), U)          
      - fvc::div(nuTensor() & dev2(T(fvc::grad(U))), "div((nuEff*dev(T(grad(U)))))")
    );
}


tmp<fvVectorMatrix> laminarTwoFluid::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{        
    // Viscosity tensor
    const volSymmTensorField muTensor("muTensor", rho*nuTensor());
    
    return
    (    
      - fvm::laplacian(muTensor, U)
      - fvc::div(muTensor & dev2(T(fvc::grad(U))), "div((muEff*dev(T(grad(U)))))")
    );
}


void laminarTwoFluid::correct()
{
    if (isA<incompressibleTwoFluidMixture>(transportModel_))
    {
        incompressibleTwoFluidMixture& itfm 
            = refCast<incompressibleTwoFluidMixture>(transportModel_);
            
        // Cell gradient of alpha
        volVectorField gradAlpha(fvc::grad(itfm.alpha1()));
        
        // Face unit interface normal
        volVectorField nHat(gradAlpha/(mag(gradAlpha)+ deltaN_));

        // Interfacial viscosity tensor        
        nuTensor_ = (symmTensor::I-sqr(nHat)) * itfm.nu() + sqr(nHat) * itfm.nuInv();
    }
    else
    {
        // Interfacial viscosity tensor
        nuTensor_ = (symmTensor::I) * transportModel_.nu();
    }    
    turbulenceModel::correct();
}


bool laminarTwoFluid::read()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
