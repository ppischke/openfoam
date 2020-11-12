{
    const word alphaScheme("div(phi,alpha)");    

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    
    // Face flux compression velocity
    surfaceScalarField phic
    (
        4.0 * mag(phi)
    );
                      
    forAll(phic.boundaryField(), patchi)
    {
        fvsPatchScalarField& phicp = phic.boundaryField()[patchi];

        if (!phicp.coupled())
        {
            phicp == 0;
        }
    }
        
    // Uncompressed alpha flux
    surfaceScalarField phiAlpha
    (
        upwind<scalar>
        (
            mesh,
            phi
        ).flux(alpha1.oldTime())
    );
    
    // Solve explicitly WITHOUT MULES limiter
    MULES::explicitSolve
    (            
        geometricOneField(),
        alpha1,
        phiAlpha,
        zeroField(),
        zeroField()
    );
    alpha2 = 1.0 - alpha1;  
    
    interface.correct();
    
    // Artificial compression                
    const surfaceScalarField phir
    (
        phic * (interface.nHatf()) * mag(interface.nHatf()) * interface.cAlpha()
    );

    // Normalized gradient magnitude
    volScalarField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("V",dimVolume,0.0)
    );
    V.internalField() = mesh.V();
    V.correctBoundaryConditions();
    
    volSymmTensorField gradAlphaT
    (
        4.0 * V * inv(fvc::surfaceSum(sqr(mesh.Sf()) / mesh.magSf()))
    );   
    
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {                      
        volScalarField magGradAlpha
        (
            mag(fvc::grad(alpha1) & gradAlphaT) * (1.0 - interface.contactBoundary())
        );
        
        volScalarField maxGradAlpha
        (
            magGradAlpha
        );
                
        forAll(nei, facei)
        {
            maxGradAlpha[own[facei]] = max(maxGradAlpha[own[facei]], magGradAlpha[nei[facei]]);
            maxGradAlpha[nei[facei]] = max(maxGradAlpha[nei[facei]], magGradAlpha[own[facei]]);
        }

        // Gradient quality
        const surfaceScalarField gradAlphar
        (
            localMax<scalar>
            (
                mesh
            ).interpolate(maxGradAlpha)
        );
           
        // Compress alpha flux
        phiAlpha += upwind<scalar>
        (
            mesh,
          - upwind<scalar>
            (
                mesh,
              - phir * (1.0 - gradAlphar)
            ).flux(alpha2)
        ).flux(alpha1);
        			
        // Solve explicitly WITHOUT MULES limiter
        MULES::explicitSolve
        (            
            geometricOneField(),
            alpha1,
            phiAlpha,
            zeroField(),
            zeroField()
        );
        alpha2 = 1.0 - alpha1;  
    }
    interface.correct();

    rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
