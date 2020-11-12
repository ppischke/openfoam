/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.3.x                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0.000000e+00";
    object      alpha.liq;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    bottom
    {
        type            fixedValue;
        value           uniform 1;
    }
    left
    {
        type            constantAlphaContactAngle;
        limit           alpha;
        theta0          45;
        value           uniform 0;
    }
    right
    {
        type            constantAlphaContactAngle;
        limit           alpha;
        theta0          45;
        value           uniform 0;
    }
    top
    {
        type            fixedValue;
        value           uniform 0;
    }
    backFront
    {
        type            empty;
    }
}


// ************************************************************************* //
