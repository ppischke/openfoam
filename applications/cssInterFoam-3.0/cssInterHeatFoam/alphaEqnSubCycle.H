if (nAlphaSubCycles > 1)
{
    const dimensionedScalar totalDeltaT = runTime.deltaT();
    
    surfaceScalarField rhoPhiSum
    (
        IOobject
        (
            "rhoPhiSum",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("0", rhoPhi.dimensions(), 0)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;

// ec: calculating thermal properties

//kEff==1/(Pr1*alpha1/Cp1/rho1/nu1+Pr2*alpha2/Cp2/rho2/nu2);
kEff==nuInv/(Pr1*alpha1/Cp1/rho1+Pr2*alpha2/Cp2/rho2);

// end ec
