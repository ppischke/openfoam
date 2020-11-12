        // Solve energy equation
 //       fvScalarMatrix TEqn
 //       (
 //           fvm::ddt(rho*Cp, T)
 //         + fvm::div(fvc::interpolate(Cp)*rhoPhi, T)
 //         - fvm::laplacian(k, T)
 //       );

        

  	fvScalarMatrix TEqn
        (
            fvm::ddt((Cp1*rho1*alpha1+Cp2*rho2*alpha2), T)
          + fvm::div(phi*fvc::interpolate(Cp1*rho1*alpha1+Cp2*rho2*alpha2), T)
          - fvm::laplacian(kEff, T)
        );

        TEqn.solve();
	
	
