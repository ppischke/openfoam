Info << "Write Sample" << endl;

// volume

scalar volume = sum(alpha1.internalField()*mesh.V());

// cell coordinates

volScalarField ccx = mesh.C().component(0);
volScalarField ccy = mesh.C().component(1);
volScalarField ccz = mesh.C().component(2);

// mass centre

scalar cx = sum(alpha1.internalField()*ccx.internalField()*mesh.V())/sum(alpha1.internalField()*mesh.V());
scalar cy = 0; // symmetry plane
scalar cz = sum(alpha1.internalField()*ccz.internalField()*mesh.V())/sum(alpha1.internalField()*mesh.V());

// second moment

scalar c2x = sum(alpha1.internalField()*pow((ccx.internalField()-cx),2)*mesh.V())/sum(alpha1.internalField()*mesh.V());
scalar c2y = sum(alpha1.internalField()*pow((ccy.internalField()-cy),2)*mesh.V())/sum(alpha1.internalField()*mesh.V());
scalar c2z = sum(alpha1.internalField()*pow((ccz.internalField()-cz),2)*mesh.V())/sum(alpha1.internalField()*mesh.V());

// mean velocities

scalar xVel = sum(alpha1.internalField()*U.component(0)*rho1.value()*mesh.V())/(rho1.value()*volume);
scalar yVel = sum(alpha1.internalField()*U.component(1)*rho1.value()*mesh.V())/(rho1.value()*volume);
scalar zVel = sum(alpha1.internalField()*U.component(2)*rho1.value()*mesh.V())/(rho1.value()*volume);

// mean pressure

scalar pMeanIn = sum(pos(alpha1.internalField()-0.999)*p_rgh.internalField()*mesh.V())/(1e-20+sum(pos(alpha1.internalField()-0.999)*mesh.V()));

scalar pMeanOut = sum(pos(-alpha1.internalField()+0.001)*p_rgh.internalField()*mesh.V())/(1e-20+sum(pos(-alpha1.internalField())*mesh.V()));

// contact patchi

// contact area, mean pressure

scalar contArea = 0;
scalar counter=0;

forAll(mesh.boundary(), patchi)
    {
        if (isA<alphaContactAngleFvPatchScalarField>(alpha1.boundaryField()[patchi]))
        {
			// area
			
			const scalarField& Ap = mesh.magSf().boundaryField()[patchi];
			
			const scalarField& alphaBoundary = mag(alpha1.boundaryField()[patchi]);
			
			contArea += sum(Ap*alphaBoundary);
			
		}
		
	}

// axes

// cell list axis 

DynamicList<scalar> alphaXPosweight;
DynamicList<scalar> cellXPos;
DynamicList<scalar> alphaXNegweight;
DynamicList<scalar> cellXNeg;
DynamicList<scalar> alphaYPosweight;
DynamicList<scalar> cellYPos;
// DynamicList<scalar> alphaYNegweight;  // symmetry
// DynamicList<scalar> cellYNeg;
DynamicList<scalar> alphaZPosweight;
DynamicList<scalar> cellZPos;
DynamicList<scalar> alphaZNegweight;
DynamicList<scalar> cellZNeg;

scalar xPos;
scalar xNeg;
scalar yPos;
scalar yNeg;
scalar zPos;
scalar zNeg;

scalar cellHeight;

int counterXpos=0;
int counterXneg=0;
int counterYpos=0;
// int counterYneg=0;
int counterZpos=0;
int counterZneg=0;

const faceList & ff = mesh.faces(); 
const pointField & pp = mesh.points();

// loop to find cells on axis

forAll(mesh.C(), celli)
	{
		const cell & cc = mesh.cells()[celli];
		labelList pLabels(cc.labels(ff));
  		pointField pLocal(pLabels.size(), vector::zero);
	
		forAll (pLabels, pointi)
			pLocal[pointi] = pp[pLabels[pointi]];

			//Determining the cells in x pos direction
			if(Foam::min(pLocal & vector(0,1,0))<=cy && Foam::max(pLocal & vector(0,1,0))>=cy 
			&& Foam::min(pLocal & vector(0,0,1))<=cz && Foam::max(pLocal & vector(0,0,1))>=cz)
			{
				
				if (mesh.C()[celli].component(0) > cx)
				{
				alphaXPosweight(counterXpos)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				cellHeight = Foam::max(pLocal & vector(1,0,0))-Foam::min(pLocal & vector(1,0,0));
				cellXPos(counterXpos)= mesh.C()[celli].component(0)+((alpha1[celli])-0.5)*cellHeight;
				counterXpos++;		
				}
				else
				{
				alphaXNegweight(counterXneg)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				cellHeight = Foam::max(pLocal & vector(1,0,0))-Foam::min(pLocal & vector(1,0,0));
				cellXNeg(counterXneg)= mesh.C()[celli].component(0)-((alpha1[celli])-0.5)*cellHeight;
				counterXneg++;
				}

			}	
		
			//Determining the cells in y direction
			if(Foam::min(pLocal & vector(1,0,0))<=cx && Foam::max(pLocal & vector(1,0,0))>=cx 
			&& Foam::min(pLocal & vector(0,0,1))<=cz && Foam::max(pLocal & vector(0,0,1))>=cz)
			{
				if (mesh.C()[celli].component(1) > cy)
				{
				alphaYPosweight(counterYpos)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				cellHeight = Foam::max(pLocal & vector(0,1,0))-Foam::min(pLocal & vector(0,1,0));
				cellYPos(counterYpos)= mesh.C()[celli].component(1)+((alpha1[celli])-0.5)*cellHeight;
				counterYpos++;
				}
				// else
				// {
				// alphaYNegweight(counterYneg)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				// cellHeight = Foam::max(pLocal & vector(0,1,0))-Foam::min(pLocal & vector(0,1,0));
				// cellYNeg(counterYneg)= mesh.C()[celli].component(1)-((alpha1[celli])-0.5)*cellHeight;
				// counterYneg++;
				// }

			}
			
			//Determining the cells in z direction
			if(Foam::min(pLocal & vector(0,1,0))<=cy && Foam::max(pLocal & vector(0,1,0))>=cy 
			&& Foam::min(pLocal & vector(1,0,0))<=cx && Foam::max(pLocal & vector(1,0,0))>=cx)
			{
				if (mesh.C()[celli].component(2) > cz)
				{
				alphaZPosweight(counterZpos)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				cellHeight = Foam::max(pLocal & vector(0,0,1))-Foam::min(pLocal & vector(0,0,1));
				cellZPos(counterZpos)= mesh.C()[celli].component(2)+((alpha1[celli])-0.5)*cellHeight;
				counterZpos++;
				}
				else
				{
				alphaZNegweight(counterZneg)=mag(2*(0.5-mag(alpha1[celli]-0.5)));
				cellHeight = Foam::max(pLocal & vector(0,0,1))-Foam::min(pLocal & vector(0,0,1));
				cellZNeg(counterZneg)= mesh.C()[celli].component(2)-((alpha1[celli])-0.5)*cellHeight;
				counterZneg++;
				}

			}
	}
			
    // find 3 interface cells, interpolate surface position
	
	labelList order; // track index after sorting 
	sortedOrder(alphaXPosweight, order);
	sort(alphaXPosweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	label maxIndexOrder=order.size()-1;  
	
	
	// x axis
	
	if (max(alphaXPosweight) < 0.001)  // liquid is on wall
	{
		xPos = max(cellXPos)-cx;
	}
	else
	{
	xPos = (alphaXPosweight[maxIndexOrder]*cellXPos[order[maxIndexOrder]]
			+alphaXPosweight[maxIndexOrder-1]*cellXPos[order[maxIndexOrder-1]]
			+alphaXPosweight[maxIndexOrder-2]*cellXPos[order[maxIndexOrder-2]]
			+alphaXPosweight[maxIndexOrder-3]*cellXPos[order[maxIndexOrder-3]]
			+alphaXPosweight[maxIndexOrder-4]*cellXPos[order[maxIndexOrder-4]])
			/(alphaXPosweight[maxIndexOrder]+alphaXPosweight[maxIndexOrder-1]+alphaXPosweight[maxIndexOrder-2]+alphaXPosweight[maxIndexOrder-3]+alphaXPosweight[maxIndexOrder-4]+1e-20)-cx;
	};
	
	sortedOrder(alphaXNegweight, order);
	sort(alphaXNegweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	maxIndexOrder=order.size()-1;
	

	if (max(alphaXNegweight) < 0.001)
	{
		xNeg = -min(cellXNeg)+cx;
	}
	else
	{
	xNeg = -(alphaXNegweight[maxIndexOrder]*cellXNeg[order[maxIndexOrder]]
		+alphaXNegweight[maxIndexOrder-1]*cellXNeg[order[maxIndexOrder-1]]
		+alphaXNegweight[maxIndexOrder-2]*cellXNeg[order[maxIndexOrder-2]]
		+alphaXNegweight[maxIndexOrder-3]*cellXNeg[order[maxIndexOrder-3]]
		+alphaXNegweight[maxIndexOrder-4]*cellXNeg[order[maxIndexOrder-4]])
		/(alphaXNegweight[maxIndexOrder]+alphaXNegweight[maxIndexOrder-1]+alphaXNegweight[maxIndexOrder-2]+alphaXNegweight[maxIndexOrder-3]+alphaXNegweight[maxIndexOrder-4]+1e-20)+cx;
	};
	
	
	// y axis
	
	sortedOrder(alphaYPosweight, order);
	sort(alphaYPosweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	maxIndexOrder=order.size()-1;  
	
	
	if (max(alphaYPosweight) < 0.001)
	{
		yPos = max(cellYPos)-cy;
	}
	else
	{
	yPos = (alphaYPosweight[maxIndexOrder]*cellYPos[order[maxIndexOrder]]
			+alphaYPosweight[maxIndexOrder-1]*cellYPos[order[maxIndexOrder-1]]
			+alphaYPosweight[maxIndexOrder-2]*cellYPos[order[maxIndexOrder-2]]
			+alphaYPosweight[maxIndexOrder-3]*cellYPos[order[maxIndexOrder-3]]
			+alphaYPosweight[maxIndexOrder-4]*cellYPos[order[maxIndexOrder-4]])
			/(alphaYPosweight[maxIndexOrder]+alphaYPosweight[maxIndexOrder-1]+alphaYPosweight[maxIndexOrder-2]+alphaYPosweight[maxIndexOrder-3]+alphaYPosweight[maxIndexOrder-4]+1e-20)-cy;
	};
	
	// sortedOrder(alphaYNegweight, order);
	// sort(alphaYNegweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	// maxIndexOrder=order.size()-1;
	

	// if (max(alphaYNegweight) < 0.001)
	// {
		// yNeg = -min(cellXNeg)+cy;
	// }
	// else
	// {
	// yNeg = -(alphaYNegweight[maxIndexOrder]*cellYNeg[order[maxIndexOrder]]
		// +alphaYNegweight[maxIndexOrder-1]*cellYNeg[order[maxIndexOrder-1]]
		// +alphaYNegweight[maxIndexOrder-2]*cellYNeg[order[maxIndexOrder-2]]
		// +alphaYNegweight[maxIndexOrder-3]*cellYNeg[order[maxIndexOrder-3]]
		// +alphaYNegweight[maxIndexOrder-4]*cellYNeg[order[maxIndexOrder-4]])
		// /(alphaYNegweight[maxIndexOrder]+alphaYNegweight[maxIndexOrder-1]+alphaYNegweight[maxIndexOrder-2]+alphaYNegweight[maxIndexOrder-3]+alphaYNegweight[maxIndexOrder-4]+1e-20)+cy;
		
	// };
	yNeg = -yPos;
	
	// z axis
	
	sortedOrder(alphaZPosweight, order);
	sort(alphaZPosweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	maxIndexOrder=order.size()-1;  
	

	if (max(alphaZPosweight) < 0.001)
	{
		zPos = max(cellZPos)-cz;
	}
	else
	{	
	zPos = (alphaZPosweight[maxIndexOrder]*cellZPos[order[maxIndexOrder]]
			+alphaZPosweight[maxIndexOrder-1]*cellZPos[order[maxIndexOrder-1]]
			+alphaZPosweight[maxIndexOrder-2]*cellZPos[order[maxIndexOrder-2]]
			+alphaZPosweight[maxIndexOrder-3]*cellZPos[order[maxIndexOrder-3]]
			+alphaZPosweight[maxIndexOrder-4]*cellZPos[order[maxIndexOrder-4]])
			/(alphaZPosweight[maxIndexOrder]+alphaZPosweight[maxIndexOrder-1]+alphaZPosweight[maxIndexOrder-2]+alphaZPosweight[maxIndexOrder-3]+alphaZPosweight[maxIndexOrder-4]+1e-20)-cz;
			
	sortedOrder(alphaZNegweight, order);
	sort(alphaZNegweight); // max values are interface cells, use max, max-1, max-2 cells for interpolation
	maxIndexOrder=order.size()-1;
	};

	if (max(alphaZNegweight) < 0.001)
	{
		zNeg = -min(cellZNeg)+cz;
	}
	else
	{
	zNeg = -(alphaZNegweight[maxIndexOrder]*cellZNeg[order[maxIndexOrder]]
		+alphaZNegweight[maxIndexOrder-1]*cellZNeg[order[maxIndexOrder-1]]
		+alphaZNegweight[maxIndexOrder-2]*cellZNeg[order[maxIndexOrder-2]]
		+alphaZNegweight[maxIndexOrder-3]*cellZNeg[order[maxIndexOrder-3]]
		+alphaZNegweight[maxIndexOrder-4]*cellZNeg[order[maxIndexOrder-4]])
		/(alphaZNegweight[maxIndexOrder]+alphaZNegweight[maxIndexOrder-1]+alphaZNegweight[maxIndexOrder-2]+alphaZNegweight[maxIndexOrder-3]+alphaZNegweight[maxIndexOrder-4]+1e-20)+cz;		
	};

// write Data

ofstream myfile ("SampleLog", ios_base::out | ios_base::app);

myfile << "\n";				 												
myfile << runTime.timeName() << "\t";
myfile << volume << "\t";
myfile << cx << "\t";
myfile << cy << "\t";
myfile << cz << "\t";
myfile << c2x << "\t";
myfile << c2y << "\t";
myfile << c2z << "\t";
myfile << xPos << "\t";
myfile << xNeg << "\t";
myfile << yPos << "\t";
myfile << yNeg << "\t";
myfile << zPos << "\t";
myfile << zNeg << "\t";
myfile << xVel << "\t";
myfile << yVel << "\t";
myfile << zVel << "\t";
myfile << contArea << "\t";
myfile << pMeanIn << "\t";
myfile << pMeanOut << "\t";

myfile.close();
