/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "partialSlipFvPatchFields.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "localMax.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::constant::mathematical::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::GeometricBoundaryField& nHatb,
    surfaceVectorField::GeometricBoundaryField& gradAlphab,
    surfaceScalarField::GeometricBoundaryField& snGradAlphab
)
{
    const fvMesh& mesh = alpha1_.mesh();
    const volScalarField::GeometricBoundaryField& abf = alpha1_.boundaryField();
    const volVectorField::GeometricBoundaryField& ubf = U_.boundaryField();
    
    const fvBoundaryMesh& boundary = mesh.boundary();

    // Unmark contact boundary cells
    contactBoundary_ = 0.0;

    forAll(boundary, patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        {
            alphaContactAngleFvPatchScalarField& acap =
                const_cast<alphaContactAngleFvPatchScalarField&>
                (
                    refCast<const alphaContactAngleFvPatchScalarField>
                    (
                        abf[patchi]
                    )
                );
                
            fvsPatchVectorField& nHatf = nHatb[patchi];
            fvsPatchVectorField& gradAlphaf = gradAlphab[patchi];                                    
            fvsPatchScalarField& snGradAlpha = snGradAlphab[patchi];
                                    
            // Contact angle
            const scalarField theta = convertToRad*acap.theta(U_.boundaryField()[patchi], nHatf);

            // Boundary face normal
            const vectorField nf = boundary[patchi].nf();
            
            // Split gradient of alpha
            vectorField gradAlphat = gradAlphaf - nf * (nf & gradAlphaf);
            vectorField gradAlphan = nf * mag(gradAlphat) / tan(theta);
                                                    
            // Face gradient of alpha
            gradAlphaf = gradAlphan + gradAlphat;
                        
            // Face unit interface normal
            nHatf = gradAlphaf/(mag(gradAlphaf)+deltaN_.value()); 

            // Surface normal gradient
            snGradAlpha = 3.0/2.0*gradAlphaf & nf;
            
            acap.gradient() = snGradAlpha;
            acap.evaluate();            

            // Mark contact boundary cells
            const unallocLabelList& contactCells =
                boundary[patchi].faceCells();
      
            forAll(boundary[patchi], facei)
            {
                contactBoundary_[contactCells[facei]] = 1.0;
                contactBoundary_.boundaryField()[patchi][facei] = 1.0;
            }            
        }
        
        
        if (isA<partialSlipFvPatchVectorField>(ubf[patchi]))
        {
            partialSlipFvPatchVectorField& psp =
                const_cast<partialSlipFvPatchVectorField&>
                (
                    refCast<const partialSlipFvPatchVectorField>
                    (
                        ubf[patchi]
                    )
                );
             
            const fvPatchScalarField& alphaf = abf[patchi];
             
            psp.valueFraction() = scalar(1) - scalar(4) * (alphaf) * (scalar(1)-alphaf);  
            psp.evaluate();
        }
    }
}


void Foam::interfaceProperties::calculateK()
{
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField nf(Sf/mag(Sf));

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "alpha"));        
    
    // Cell unit interface normal
    const volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
        
    // Face gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));    
    
    // Face unit interface normal    
    surfaceVectorField nHatf(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Surface-normal gradient of alpha
    surfaceScalarField snGradAlpha(fvc::snGrad(alpha1_, "alpha"));

    // Correct boundary conditions
    correctContactAngle(nHatf.boundaryField(), gradAlphaf.boundaryField(), snGradAlpha.boundaryField());

    // Face unit interface normal flux
    nHatf_ = nHatf & nf;
     
    // Interface tangential expression for curvature
    Kf_ = tr(fvc::interpolate(fvc::grad(nHatf)) & (sqr(nHatf) - symmTensor::I));
      
    // Continuous surface force approach
    phii_ = sigma_ * Kf_ * snGradAlpha;
}


void Foam::interfaceProperties::calculateSurfaceStress()
{    
    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const surfaceVectorField nf(Sf/mag(Sf));

    // Contact boundary
    const surfaceScalarField contactBoundaryf
    (
        localMax<scalar>(mesh).interpolate(contactBoundary_)
    );

    // Cell gradient of alpha
    const volVectorField gradAlpha(fvc::grad(alpha1_, "alpha"));        
    
    // Cell unit interface normal
    const volVectorField nHat(gradAlpha/(mag(gradAlpha) + deltaN_));
        
    // Face gradient of alpha
    surfaceVectorField gradAlphaf(fvc::interpolate(gradAlpha));    
    
    // Face unit interface normal    
    surfaceVectorField nHatf(gradAlphaf/(mag(gradAlphaf) + deltaN_));

    // Surface-normal gradient of alpha
    surfaceScalarField snGradAlpha(fvc::snGrad(alpha1_, "alpha"));

    // Correct boundary conditions
    correctContactAngle(nHatf.boundaryField(), gradAlphaf.boundaryField(), snGradAlpha.boundaryField());
        
    // Face unit interface normal flux
    nHatf_ = nHatf & nf;        
    
    // Interface tangential expression for curvature
    Kf_ = tr(fvc::interpolate(fvc::grad(nHatf)) & (sqr(nHatf) - symmTensor::I));

    // Interface surface stress        
    sStress_ = sigmaf_ * mag(gradAlphaf) * (symmTensor::I - sqr(nHatf)) & mesh.Sf();    
        
    // Continuous surface stress approach, blend into surface force at boundary
    phii_ = (1.0-contactBoundaryf) * (fvc::interpolate(fvc::div(sStress_)) & nf)
              + (contactBoundaryf) * (sigmaf_ * Kf_ * snGradAlpha);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
    transportPropertiesDict_(dict),
    alphaControls_(alpha1.mesh().solverDict(alpha1.name())),
    
    sigma_
    (
        transportPropertiesDict_.lookup("sigma")
    ),

    cAlpha_
    (
        alphaControls_.lookupOrDefault<scalar>("cAlpha", 1.0)
    ),
    
    deltaN_
    (
        alphaControls_.lookupOrDefault<scalar>("deltaN", SMALL) / pow(average(alpha1.mesh().V()), 1.0/3.0)        
    ),

    alpha1_(alpha1),
    U_(U),
    

    sigmaf_
    (
        IOobject
        (
            "sigmaf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        sigma_
    ),

    nHatf_
    (
        IOobject
        (
            "nHatf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("nHatf", dimless, 0.0)
    ),

    Kf_
    (
        IOobject
        (
            "Kf",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("Kf", dimless/dimLength, 0.0)
    ),
    
    sStress_    
    (
        IOobject
        (
            "sStress",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedVector("sStress", dimensionSet(1,1,-2,0,0), vector::zero)
    ),
    
    phii_    
    (
        IOobject
        (
            "phii",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        dimensionedScalar("phii", dimensionSet(1,-2,-2,0,0), 0.0)
    ), 

    contactBoundary_
    (
        IOobject
        (
            "contactBoundary",
            alpha1_.time().timeName(),
            alpha1_.mesh()
        ),
        alpha1_.mesh(),
        scalar(0.0)
    ) 
{
    calculateSurfaceStress();
}


// ************************************************************************* //
