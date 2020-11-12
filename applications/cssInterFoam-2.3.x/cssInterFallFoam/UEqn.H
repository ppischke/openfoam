    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)   
      + turbulence->divDevRhoReff(rho, U)   
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==            
            fvc::reconstruct
            (
                (
                    interface.phii()
                  + fvc::interpolate(rho-rho2)*gf
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );
    }
