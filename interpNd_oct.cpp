#include <octave/oct.h>
//#include <octave/ov.h>
//#include <dMatrix.h>

/* C++ includes */
#include <math.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <cmath>
using namespace std;

/* interpolation includes */
#include "src/interpND.h"


void interpolate_main(double* v0, double p[], double* X, double* Y, double* Z, double* V,
	int Nx, int Ny, int Nz, int dim, int Np, bool LayElev, int no_data, int imeth){
	double **XX, ***VV, ***ZZ;

	int matrix[4]; matrix[0] = dim; matrix[1] = Nz;
	matrix[2] = Ny; matrix [3] = Nx;

	/* allocate space for matrices */
    XX = new double*[dim];
    for (int i = 0; i < dim; i++)
        XX[i] = new double [matrix[3-i]];

    VV = new double**[matrix[1]];
    for (int i = 0; i < matrix[1]; i++){
        VV[i] = new double*[matrix[2]];
        for (int j = 0; j < matrix[2]; j++)
            VV[i][j] = new double[matrix[3]];  
    }
    if (LayElev){
    	ZZ = new double**[matrix[1]];
	    for (int i = 0; i < matrix[1]; i++){
	        ZZ[i] = new double*[matrix[2]];
	        for (int j = 0; j < matrix[2]; j++)
	            ZZ[i][j] = new double[matrix[3]];  
	    }
    }

    // put values into interpolant matrices
    int cnt = 0;
    for (int k = 0; k < matrix[1]; k++)
    	for (int j = 0; j < matrix[3]; j++)
    		for (int i = 0; i < matrix[2]; i++){
    			VV[k][i][j] = V[cnt];
    			//cout << VV[k][i][j] << endl;
    			cnt++;
    		}
    if (dim >= 1){
        for (int i = 0; i < matrix[3]; i++){
            XX[0][i] = X[i];
            //cout << XX[0][i] << endl;
        }
    }
    if (dim >= 2){
        for (int i = 0; i < matrix[2]; i++){
            XX[1][i] = Y[i];
            //cout << XX[1][i] << endl;
        }
    }
    if (dim >= 3){
    	if (LayElev == false){
        	for (int i = 0; i < matrix[1]; i++){
            	XX[2][i] = Z[i];
            	//cout << XX[2][i] << endl;
        	}
        }
        else{
        	for (int k = 0; k < matrix[1]; k++)
		    	for (int j = 0; j < matrix[3]; j++)
		    		for (int i = 0; i < matrix[2]; i++){
		    			ZZ[k][i][j] = Z[cnt];
		    			cout << ZZ[k][i][j] << endl;
		    			cnt++;
		    		}

        }
    }
    


    double *x0 = new double[dim];



	if (dim == 1){
		interp<1> myinterp;
		myinterp.no_data = no_data;
		myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
		myinterp.set_Method(imeth);
    	for (int i = 0; i < Np; i++){
			x0[0] = p[i];
			v0[i] = myinterp.interpolate(x0);
		}
	}
	else if (dim == 2){
		interp<2> myinterp;
		myinterp.no_data = no_data;
		myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
		myinterp.set_Method(imeth);
		for (int i = 0; i < Np; i++){
			x0[0] = p[i];
			x0[1] = p[i + Np];
			v0[i] = myinterp.interpolate(x0);
		}

	}
	else if (dim == 3){
		interp<3> myinterp;
		myinterp.no_data = no_data;
		myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
		myinterp.set_Method(imeth);
		if (LayElev){
			myinterp.set_LayElev(true);
			myinterp.set_Z(ZZ);
		}
		for (int i = 0; i < Np; i++){
            x0[0] = p[i];
            x0[1] = p[i + Np];
            x0[2] = p[i + 2*Np];
            v0[i]=myinterp.interpolate(x0);
        }

	}
	


	// Destroy matrices
    for (int i = 0; i < dim; i++)
        delete[] XX[i];
    delete[] XX;
    
    for (int i = 0; i < matrix[1]; i++){
        for (int j = 0; j < matrix[2]; j++)
            delete[] VV[i][j];
        delete[] VV[i];
    }
    delete[] VV;
    if (LayElev){
    	for (int i = 0; i < matrix[1]; i++){
	        for (int j = 0; j < matrix[2]; j++)
	            delete[] ZZ[i][j];
	        delete[] ZZ[i];
	    }
    	delete[] ZZ;
    }
    
    delete[] x0;

}





DEFUN_DLD (interpNd_oct, args, nargout , "v0=interpNd_oct(p,X,Y,Z,V,nodata) \n \n"){


	octave_value_list retval;
	retval.resize (1);
	int nargin = args.length ();
	if (nargin != 6){
		std::cout << "Wrong number of inputs" << std::endl;
		print_usage ();
		retval(0)=0;
		return retval;
	}

	Matrix p;  //  List of points to be interpolated
	Matrix X;  // coordinates along X direction
	Matrix Y;  // coordinates along Y direction
	NDArray Z;  // coordinates along Z direction
	NDArray V;  // interpolation data
	Matrix Nodata; // nodata value
	Matrix method; // method 1-> linear, 2-> nearest

	//assign input arguments to matrices
	p  =  args(0).matrix_value();
	X  =  args(1).matrix_value();
	Y  =  args(2).matrix_value();
	Z  =  args(3).array_value();
	V  =  args(4).array_value();
	Nodata  =  args(5).matrix_value();
	int no_data = (int)Nodata(0);
	method = args(6).matrix_value();
	int imeth = (int)method(0);

	/* find the number of points to interpolate */
	dim_vector size_Dim;
	size_Dim = p.dims();
	int Np = size_Dim(0);
	int dim = size_Dim(1);

	size_Dim = X.dims();
	int Nx  = size_Dim(0);
	size_Dim = Y.dims();
	int Ny  = size_Dim(0); if (Ny == 0) Ny = 1; 
	size_Dim = Z.dims();
	int Nz  = size_Dim(0);  if (Nz == 0) Nz = 1;
	bool LayElev = false;
	if (size_Dim.length() == 3)
		LayElev = true;


	//cout << size_Dim.length() << endl;
 	cout << Np << " " << dim << endl;
 	cout << Ny << " " << Nz << endl;

 	Matrix mV0(Np,1);
	
	interpolate_main(mV0.fortran_vec(), p.fortran_vec(), X.fortran_vec(),
					 Y.fortran_vec(), Z.fortran_vec(), V.fortran_vec(),
					 Nx, Ny, Nz, dim, Np, LayElev, no_data, imeth);
	
	retval(0) = mV0;
	return octave_value_list(retval);

}