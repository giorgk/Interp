//============================================================================
// Name        : Interpolation.cpp
// Author      : George Kourakos
// Version     :
// Copyright   :
// Description : Example script on how to use the interpolation class
//============================================================================

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <cmath>
using namespace std;

#include "interpND.h"


int main() {
	cout << "@@!!!Hello World!!!@@" << endl; // prints @@!!!Hello World!!!@@


/*
	//-----------1D interpolation -----------
	char* namefile = "interp1D_example.txt";
	interp<1> myinterp1;
	myinterp1.get_data_file(namefile);
	myinterp1.set_Method(2);
	double X1[1];
	X1[0] = -1.0;
	while (X1[0] < 11){
		cout << X1[0] << "\t" << myinterp1.interpolate(X1) << endl;
		X1[0] = X1[0] + 0.33;
	}
	myinterp1.destroy();
*/


	//------------2D interpolation ------------
	char* namefile1 = "interp2D_example.txt";
	interp<2> myinterp2;
	myinterp2.set_nP(11);
	myinterp2.get_data_file(namefile1);
	myinterp2.set_Method(2);
	double X2[2];
	X2[0] = -3.5; X2[1] = -3.5;
	while (X2[1] < 3.5){
		while (X2[0] < 3.5){
			cout << X2[0] << "\t" << X2[1] << "\t" << myinterp2.interpolate(X2) << endl;
			X2[0] = X2[0] + 0.1;
		}
		X2[0] = -3.5;
		X2[1] = X2[1] + 0.1;
	}
	myinterp2.destroy();


/*
	//------------3D interpolation ------------
	char* namefile2 = "interp3D_example.txt";
	interp<3> myinterp3;
	myinterp3.get_data_file(namefile2);
	myinterp3.set_Method(2);
	double X3[3];
	X3[0] = -1.0; X3[1] = -4.5; X3[2] = -4;
	while (X3[2] < 4.5){
		while (X3[1] < 4.5){
			while (X3[0] <12){
				cout << X3[0] << "\t" << X3[1] << "\t" << X3[2] << "\t" << myinterp3.interpolate(X3) << endl;
				X3[0] = X3[0] + 1;
			}
			X3[0] = -1.0;
			X3[1] = X3[1] + 0.7;
		}
		X3[1] = -4.5;
		X3[2] = X3[2] + 0.7;
	}
	myinterp3.destroy();



	// ------------3D interpolation (layers------------
	char* namefile3 = "interp3D_example1.txt";
	interp<3> myinterp4;
	myinterp4.set_LayElev(true);
	myinterp4.get_data_file(namefile3);
	myinterp4.set_Method(2);
	double X4[3];
	X4[0] = -1.0; X4[1] = -3.5; X4[2] = -4.5;
	while (X4[2] < 6.2){
		while (X4[1] < 3.5){
			while (X4[0] < 11){
				cout << X4[0] << "\t" << X4[1] << "\t" << X4[2] << "\t" << myinterp4.interpolate(X4) << endl;
				X4[0] = X4[0] + 1;
			}
			X4[0] = -1.0;
			X4[1] = X4[1] + 0.7;
		}
		X4[1] = -3.5;
		X4[2] = X4[2] + 0.7;
	}
	myinterp4.destroy();



*/
	return 0;
}
