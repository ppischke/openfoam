/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::interfaceProperties

Description
    Contains the interface properties.

    Properties to aid interFoam:
    -# Correct the alpha boundary condition for dynamic contact angle.
    -# Calculate interface curvature.

SourceFiles
    interfaceProperties.C

\*---------------------------------------------------------------------------*/

#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

        //- Keep a reference to the transportProperties dictionary
        const dictionary& transportPropertiesDict_;
        const dictionary& alphaControls_;
        
        //- Surface tension
        const dimensionedScalar sigma_;

        //- Compression coefficient
        const scalar cAlpha_;
        
        //- Stabilisation for normalisation of the interface normal
        const dimensionedScalar deltaN_;

        const volScalarField& alpha1_;
        const volVectorField& U_;
        surfaceScalarField sigmaf_;
        surfaceScalarField nHatf_;
        surfaceScalarField Kf_;
        surfaceVectorField sStress_;
        surfaceScalarField phii_;
        volScalarField contactBoundary_;
        

    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        interfaceProperties(const interfaceProperties&);
        void operator=(const interfaceProperties&);

        //- Correction for the boundary condition on the unit normal nHat on
        //  walls to produce the correct contact dynamic angle
        //  calculated from the component of U parallel to the wall
        void correctContactAngle
        (
            surfaceVectorField::GeometricBoundaryField& nHatb,
            surfaceVectorField::GeometricBoundaryField& gradAlphab,
            surfaceScalarField::GeometricBoundaryField& snGradAlphab
        );
        
        //- Re-calculate the interface curvature
        void calculateK();
        void calculateSurfaceStress();


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions

        scalar cAlpha() const
        {
            return cAlpha_;
        }

        const dimensionedScalar& deltaN() const
        {
            return deltaN_;
        }

        const volScalarField& contactBoundary() const
        {
            return contactBoundary_;
        }
        
        const surfaceScalarField& nHatf() const
        {
            return nHatf_;
        }

        const surfaceScalarField& Kf() const
        {
            return Kf_;
        }

        const surfaceVectorField& sStress() const
        {
            return sStress_;
        }
        
        const surfaceScalarField& phii() const
        {
            return phii_;
        }        

        const dimensionedScalar& sigma() const
        {
            return sigma_;
        }
        
        surfaceScalarField& sigmaf()
        {
            return sigmaf_;
        }
        
        void correct()
        {
            calculateSurfaceStress();
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
