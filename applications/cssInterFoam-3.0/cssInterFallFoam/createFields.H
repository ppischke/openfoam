    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"


    Info<< "Reading transportProperties\n" << endl;
    incompressibleTwoFluidMixture twoPhaseProperties(U, phi);

    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());

    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        alpha1*rho1 + alpha2*rho2,
        alpha1.boundaryField().types()
    );
    rho.oldTime();


    // Mass flux
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho)*phi
    );
    
    
    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, U, twoPhaseProperties);


    // Construct incompressible turbulence model
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, twoPhaseProperties)
    );

    #include "readGravitationalAcceleration.H"

    /*
    dimensionedVector g0(g);

    // Read the data file and initialise the interpolation table
    interpolationTable<vector> timeSeriesAcceleration
    (
        runTime.path()/runTime.caseConstant()/"acceleration.dat"
    );
    */
    
    Info<< "Calculating field g\n" << endl;    
    surfaceScalarField gf("ghf", g & mesh.Sf()/mesh.magSf());

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p_rgh,
        mesh.solutionDict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    fv::IOoptionList fvOptions(mesh);
