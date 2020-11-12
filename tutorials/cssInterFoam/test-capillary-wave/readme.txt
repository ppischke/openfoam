OpenFOAM tutorial for cssInterFoam: two-dimensional capillary wave

How to:
- Run blockMesh
- Run writeCellCentres
- Run setFields.m in MATLAB
- Run cssInterFoam

Validate:
- Include debug files into cssInterFoam:
    writeAxisHeader.H
    writeEnergyHeader.H
    writeAxis.H
    writeEnergy.H
  The headers are found in the parent folder and
  currently commented from the cssInterFoam code
- Recompile cssInterFoam
- Run cssInterFoam
- Run validate.m in MATLAB
  
