/* Matlab Include */
#include "mex.h"
#include "matrix.h"


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

/* Input Arguments */
#define	P_IN	    prhs[0] // Points to interpolate
#define	X_IN	    prhs[1] // X coordinates
#define	Y_IN	    prhs[2] // Y coordinates
#define	Z_IN	    prhs[3] // Z coordinates
#define	V_IN	    prhs[4] // Values
#define	NODATA_IN	prhs[5] // nodata value


/* Output Arguments */
#define	V0_OUT plhs[0]


void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray	*prhs[]){
    double *p, *X, *Y, *Z, *V, *v0;
    bool LayElev;
    double **XX, ***VV, ***ZZ, *nodata;
    int Nx, Ny, Nz;
    int matrix[4];
            
    p = mxGetPr( P_IN );
    X = mxGetPr( X_IN );
    Y = mxGetPr( Y_IN );
    Z = mxGetPr( Z_IN );
    V = mxGetPr( V_IN );
    nodata = mxGetPr( NODATA_IN );
    
    
    /* find the number of points to interpolate */
    int Np = mxGetM(P_IN);
    /* find the dimension of interpolation */
    int dim = mxGetN(P_IN); matrix[0] = dim;
    if (dim == 1){
        Nx = mxGetM(X_IN); matrix[3] = Nx; 
        matrix[2] = 1; matrix[1] = 1;
    }
    if (dim == 2){
        Nx = mxGetM(X_IN); matrix[3] = Nx;
        Ny = mxGetM(Y_IN); matrix[2] = Ny;
        matrix[1] = 1;
    }
    if (dim >= 3){
        Nx = mxGetM(X_IN); matrix[3] = Nx;
        Ny = mxGetM(Y_IN); matrix[2] = Ny;
        Nz = mxGetM(Z_IN); matrix[1] = Nz;
    }
    
    //cout << Nx << " " << Ny << " " << Nz << endl;

    
    /* find the type of interpolation in case of 3D */
    const mwSize *dims_Z;
    if (dim == 3){
        if (mxGetNumberOfDimensions(Z_IN) == 3)
            LayElev = true;
        dims_Z = mxGetDimensions(Z_IN);
    }
    //std::cout << mxGetNumberOfDimensions(Z_IN) << std::endl;
    
    std::cout << Np << " " << dim << std::endl;
    
    
    /* Create matrices for the return arguments */
	V0_OUT = mxCreateDoubleMatrix( Np, 1, mxREAL );
    
    /*  create a C pointer to a copy of the output matrix */
    v0 = mxGetPr(V0_OUT);
    
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
            int cnt = 0;
            for (int k = 0; k < matrix[1]; k++)
                for (int j = 0; j < matrix[3]; j++)
                    for (int i = 0; i < matrix[2]; i++){
                        ZZ[k][i][j] = Z[cnt];
                        //cout << ZZ[k][i][j] << endl;
                        cnt++;
                    }
            
        }
    }
    
    
    
    /* create the interpolants */
    double *x0 = new double[dim];
    if (dim == 1){
        interp<1> myinterp;
        myinterp.no_data = nodata[0];
        cout << myinterp.no_data << endl;
        myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
        for (int i = 0; i < Np; i++){
             x0[0] = p[i];
             v0[i]=myinterp.interpolate(x0);
             //cout << v0[i] << endl;
        }
    }
    else if (dim == 2){
        interp<2> myinterp;
        myinterp.no_data = nodata[0];
        myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
        for (int i = 0; i < Np; i++){
            x0[0] = p[i];
            x0[1] = p[i + Np];
             v0[i]=myinterp.interpolate(x0);
        }
    }
    else if (dim ==3){
        interp<3> myinterp;
        myinterp.no_data = nodata[0];
        myinterp.get_data_c(matrix[3], matrix[2], matrix[1], XX, VV);
        if (LayElev){
            myinterp.set_LayElev(true);
            myinterp.set_Z(ZZ);
        }
        for (int i = 0; i < Np; i++){
            x0[0] = p[i];
            x0[1] = p[i + Np];
            x0[2] = p[i + 2*Np];
            //cout << x0[0] << " " << x0[1] << " " << x0[2] << " " <<endl;
        //cout << myinterp.interpolate(x0) << endl;
            v0[i]=myinterp.interpolate(x0);
        }
        //myinterp.destroy();
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
    //delete[] VV;
    delete[] x0;
  
}