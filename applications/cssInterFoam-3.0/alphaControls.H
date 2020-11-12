const dictionary& alphaControls = mesh.solverDict(alpha1.name());

label nAlphaCorr
(
    alphaControls.lookupOrDefault<label>("nAlphaCorr", 3)
);

label nAlphaSubCycles
(
    alphaControls.lookupOrDefault<label>("nAlphaSubCycles", 2)
);

bool alphaOuterCorrectors
(
    alphaControls.lookupOrDefault<Switch>("alphaOuterCorrectors", false)
);
