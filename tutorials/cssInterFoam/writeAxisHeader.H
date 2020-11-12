
Info << "Write Sample Header" << endl;
	
ofstream myfile ("SampleLog", ios_base::out | ios_base::app);

// write header

myfile << "\n";
myfile << "time [s]\t";
myfile << "volume\t";
myfile << "cx\t";
myfile << "cy\t";
myfile << "cz\t";
myfile << "c2x\t";
myfile << "c2y\t";
myfile << "c2z\t";
myfile << "xPos\t";
myfile << "xNeg\t";
myfile << "yPos\t";
myfile << "yNeg\t";
myfile << "zPos\t";
myfile << "zNeg\t";
myfile << "xVel\t";
myfile << "yVel\t";
myfile << "zVel\t";
myfile << "contArea\t";
myfile << "pMeanIn\t";
myfile << "pMeanOut\t";

myfile.close();

