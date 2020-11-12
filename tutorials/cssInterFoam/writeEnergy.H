Info << "Write Sample Energy" << endl;

// kineticEnergy

scalarField  kinEnergy = 0.5*rho1.value()*mag(alpha1.internalField())*mesh.V()*magSqr(U.internalField());
scalar kinEnergySum = sum(kinEnergy);

// mean delta T

// scalar deltaT = (runTime.deltaTValue() + runTime.deltaT0Value())/2;
scalar deltaT = (runTime.deltaTValue());

// dissipation

volScalarField muEff("muEff", rho*twoPhaseProperties.nu());

volVectorField dissForce =  -fvc::laplacian(muEff, U)  -fvc::div(muEff*dev(T(fvc::grad(U))));

volScalarField dissEnergy = (dissForce & U);

dissEnergySum =  dissEnergySum + sum(mag(alpha1.internalField())*dissEnergy.internalField()*mesh.V()*deltaT);

// dissipation alternativ testing

scalar dissEnergy1 = 0;

const volVectorField gradUx=fvc::grad(U.component(0));
const volVectorField gradUy=fvc::grad(U.component(1));
const volVectorField gradUz=fvc::grad(U.component(2));
const volScalarField mu = twoPhaseProperties.mu();

forAll( mesh.C(), celli)
{
dissEnergy1 = dissEnergy1
+ 2*mesh.V()[celli]*deltaT*muEff[celli]*alpha1[celli]*((pow(gradUx[celli].component(0),2)+pow(gradUy[celli].component(1),2)+pow(gradUz[celli].component(2),2))
+ 0.5*pow(gradUy[celli].component(0)+gradUx[celli].component(1),2)
+ 0.5*pow(gradUz[celli].component(1)+gradUy[celli].component(2),2)
+ 0.5*pow(gradUx[celli].component(2)+gradUz[celli].component(0),2));
}

dissEnergySum1 = dissEnergySum1 + dissEnergy1;

// surfaceEnergy

volVectorField surfForce = fvc::reconstruct(-interface.phii()*mesh.magSf());

volScalarField surfEnergy = (surfForce & U);

surfEnergySum = surfEnergySum + sum(surfEnergy.internalField()*mesh.V()*deltaT);

// alternative surface Energy

volScalarField gradAlpha = mag(fvc::grad(alpha1));
scalar totalSurfaceArea = sum(gradAlpha.internalField()*mesh.V());

scalar surfEnergySum1 = totalSurfaceArea*interface.sigma().value();


// pressure

volVectorField presForce = fvc::reconstruct((ghf*fvc::snGrad(rho)+ fvc::snGrad(p_rgh))*mesh.magSf());

volScalarField presEnergy = (presForce & U);

presEnergySum = presEnergySum + sum(presEnergy.internalField()*mesh.V()*deltaT);

// write Data

ofstream myfileEnergy (DateiNameEnergy, ios_base::out | ios_base::app);

myfileEnergy << "\n";				 												
myfileEnergy << runTime.timeName() << "\t";
myfileEnergy << dissEnergySum << "\t";
myfileEnergy << dissEnergySum1 << "\t";
myfileEnergy << surfEnergySum << "\t";
myfileEnergy << kinEnergySum << "\t";
myfileEnergy << presEnergySum << "\t";
myfileEnergy << totalSurfaceArea << "\t";
myfileEnergy << surfEnergySum1 << "\t";

// myfileEnergy.close();
