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

Class
    Foam::incompressibleTwoFluidMixture

Description
    A two-phase incompressible transportModel

SourceFiles
    incompressibleTwoFluidMixture.C

\*---------------------------------------------------------------------------*/

#ifndef incompressibleTwoFluidMixture_H
#define incompressibleTwoFluidMixture_H

#include "incompressibleTwoPhaseMixture.H"
#include "incompressible/transportModel/transportModel.H"
#include "incompressible/viscosityModels/viscosityModel/viscosityModel.H"
#include "twoPhaseMixture.H"
#include "IOdictionary.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class incompressibleTwoFluidMixture Declaration
\*---------------------------------------------------------------------------*/

class incompressibleTwoFluidMixture
:
    public incompressibleTwoPhaseMixture
{
protected:

    // Protected data

        volScalarField nuInv_;

    // Private Member Functions

        //- Calculate and return the laminar viscosity
        void calcNu();


public:

    TypeName("incompressibleTwoFluidMixture");


    // Constructors

        //- Construct from components
        incompressibleTwoFluidMixture
        (
            const volVectorField& U,
            const surfaceScalarField& phi
        );


    //- Destructor
    virtual ~incompressibleTwoFluidMixture()
    {}


    // Member Functions
        
        //- Return inverse-averaged dynamic laminar viscosity
        tmp<volScalarField> muInv() const;

        //- Return inverse-averaged dynamic laminar viscosity
        tmp<surfaceScalarField> muInvf() const;
        
        //- Return the inverse-averaged kinematic laminar viscosity
        tmp<volScalarField> nuInv() const
        {
            return nuInv_;
        }

        //- Return the inverse-averaged laminar viscosity for patch
        virtual tmp<scalarField> nuInv(const label patchi) const
        {
            return nuInv_.boundaryField()[patchi];
        }
        
        //- Return the inverse-averaged face-interpolated kinematic laminar viscosity
        tmp<surfaceScalarField> nuInvf() const;

        //- Correct the laminar viscosity
        virtual void correct()
        {
            calcNu();
        }
	
	// added by ec
	//- Return the basic kinematic viscosities of both phases
        tmp<volScalarField> nu1() const;
	tmp<volScalarField> nu2() const;
	
	// end ec
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
