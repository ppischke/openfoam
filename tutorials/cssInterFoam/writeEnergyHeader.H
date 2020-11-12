
Info << "Write Sample Header Energy" << endl;

std::ostringstream ostr;

// IOdictionary writeSampleDict
    // (
        // IOobject
        // (
            // "writeSampleDict",
            // runTime.constant(),
            // mesh,
            // IOobject::MUST_READ,
            // IOobject::NO_WRITE
        // )
    // );
	
// create output file
	
int proNumb = Pstream::myProcNo();
ostr << proNumb;

std::string e = "SampleLogEnergy" + ostr.str();
const char* DateiNameEnergy = e.c_str();

ofstream myfileEnergy (DateiNameEnergy, ios_base::out | ios_base::app);

// declare variables

scalar dissEnergySum = 0;

scalar dissEnergySum1 = 0;

scalar surfEnergySum = 0;

scalar presEnergySum = 0;

// write header

myfileEnergy << "\n";
myfileEnergy << "time [s]   \t";
myfileEnergy << "dissEnergy   \t";
myfileEnergy << "dissEnergy1   \t";
myfileEnergy << "surfEnergy   \t";
myfileEnergy << "kineEnergy   \t";
myfileEnergy << "presEnergy   \t";
myfileEnergy << "surfArea   \t";
myfileEnergy << "surfEnergy1   \t";
myfileEnergy.close();