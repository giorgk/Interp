//! Interpolation Class
/*!
 * This class is a dimension independent class for interpolation of regularly spatial distributed data.
 * It is limited to 1-2 and 3D interpolation.
 *
 * The interpolation is based on three methods. If the unknown point is surrounded by known nodes a
 * [linear](http://en.wikipedia.org/wiki/Linear_interpolation),
 * [bilinear](http://en.wikipedia.org/wiki/Bilinear_interpolation),
 * [trilinear](http://en.wikipedia.org/wiki/Trilinear_interpolation) interpolation is used.
 * If there are data gaps (e.g in 3D there are up to 7 known nodes
 *  out of 8), then a [inverse distance weighting (IDW)](http://en.wikipedia.org/wiki/Inverse_distance_weighting)
 *   interpolation is performed with power parameter equal to 1.6.
 *   In the extreme case that all values around a point are missing then the code search for a number of
 *   known points in closest possible proximity and performs IDW interpolation
 *
 *   Known issues:
 *   - When there are missing data at the boundaries then some error is expected in extrapolation due to the influence
 *    of the inner points. During extrapolation it is better to set #nP = 1,
 *    which leads to nearest neighborhood interpolation
 *     */

template <int dim>
class interp{
public:
	//! A value for missing data
	int no_data;
	interp ();
	//! Reads the data from file
	void get_data_file (char* namefile);
	//! Calculate the interpolated value
	double interpolate (double x0[]);
	//! Destroys the matrices and frees memory
	void destroy();
	//! if true then a spatially varying grid is expected
	void set_LayElev(bool value);
	//! sets the number of points #nP
	void set_nP(int value);
	//! puts the data to interpolant
	void get_data_c(int Nx, int Ny, int Nz, double **X, double ***V);
	void set_Z(double ***Z);
private:
	double **X; /**< A pointer which holds the coordinates of the grid.
				- x[0] holds the x coordinates
			    - x[1] holds the y coordinate
			    - x[2] holds the z coordinate
			    Each of the above arrays have length equal to their corresponding grid size
			    */

	double ***V; /**< A pointer to hold the property values defined on the grid #X.
	 	 	 	 The variable is accessed as follows:\n
	 	 	 	 The value <code>V[ilay][irow][icol]</code> corresponds to <code>ilay</code> layer
	 	 	 	 <code>irow</code> row and <code>icol</code> column.
	 	 	 	 */

	double ***Z; /**< A pointer to hold the elevation of the grid nodes in 3D only.
					This is very practical when the elevations of each layer varies spatially,
					which is very common in groundwater relater applications This has similar structure to #V
	 	 	 	 */

	bool LayElev;/**< A boolean value. TRUE when the elevation is defined as 3D variable. FALSE if the Z is uniform.
	 	 	 	 */
	int nP; /**< Is the number of points to consider in no-data where IDW is used.
				By default this is 2^dim
	 	 	 */
	int Matrix_dim[4];/**< An array that contains the following info:
	 	 	 	 	 -# Dimension
	 	 	 	 	 -# Grid size in z direction
	 	 	 	 	 -# Grid size in y direction
	 	 	 	 	 -# Grid size in x direction
	 	 	 	 	 */
	//! returns the index and coordinates of the cell that the interpolated point is located
	void get_index(int I[], double x[], double X0[], int idim);
	//! returns the index and the coordinates of the Z cell given the I and J cell indices.
	void get_index_Z(int I[], double y[], int J[], double x[], int K[], double *z, double x0[]);
	//! checks how many missing values exist in the current cell
	int CheckNodataValues(double *Q, int nQ);
	//! identifies the neighbors of the I,J,K cell
	int ListNeighbors(int **&LIST);
	//! returns the appropriate filter for each dimension which is used during searching for neighbors
	int Filter(int **&FLT, int mode);
	//! test if the neighbors are valid
	bool testneighbor(int *Itest, int **LIST);
};


template <int dim>
interp<dim>::interp(){
	if (dim == 1){
		Matrix_dim[0] = dim;
		Matrix_dim[1] = 1;
		Matrix_dim[2] = 1;
	}
	else if (dim == 2){
		Matrix_dim[0] = dim;
		Matrix_dim[1] = 1;
	}
	else
		Matrix_dim[0] = dim;
	LayElev = false;
	nP = pow(2,dim);
}

template <int dim>
void interp<dim>::get_data_c(int Nx, int Ny, int Nz, double **tX, double ***tV){
	Matrix_dim[0] = dim;
	Matrix_dim[1] = Nz;
	Matrix_dim[2] = Ny;
	Matrix_dim[3] = Nx;

	X = &tX[0];
	V = &tV[0];

}

template <int dim>
void interp<dim>::get_data_file(char* namefile){
	cout << namefile << endl;
	ifstream* datafile = new ifstream(namefile);

	if (!datafile->good()){
		cout << "Can't open " << namefile << endl;
		return;
	}
	else{
		char buffer[1024];
		datafile->getline(buffer,1024);
		//X = new double ***[dim];
		//V = new double ***[dim];

		{// read grid dimensions and no data value
			istringstream inp(buffer);
			for (int i = 0; i < dim; i++){
				inp >> Matrix_dim[3 - i];
				cout << "dim" << i + 1 << "=" << Matrix_dim[3 - i] << " ";
			} cout << endl;
			inp >> no_data;
		}

		// allocate space for variables
		X = new double*[dim];
		for (int idim = 0; idim < dim; idim++)
			X[idim] = new double [Matrix_dim[3 - idim]];


		V = new double** [Matrix_dim[1]];
		if (LayElev == true) Z = new double** [Matrix_dim[1]];
		for (int ilay = 0; ilay < Matrix_dim[1]; ilay++){
			V[ilay] = new double* [Matrix_dim[2]];
			if (LayElev == true) Z[ilay] = new double* [Matrix_dim[2]];
			for (int irow = 0; irow < Matrix_dim[2]; irow++ ){
				V[ilay][irow] = new double [Matrix_dim[3]];
				if (LayElev == true) Z[ilay][irow] = new double [Matrix_dim[3]];
			}
		}


		{// read grid coordinates
			for (int idim = 0; idim < dim; idim++){
				if ((idim == 2) & (LayElev == true)){
					// Z elevation in 3D interpolation has a different format. Each value has each own elevation
					for (int ilay = 0; ilay < Matrix_dim[1]; ilay++){
						for (int irow = 0; irow < Matrix_dim[2]; irow++){
							datafile->getline(buffer,1024);
							istringstream inp(buffer);
							for (int icol = 0; icol < Matrix_dim[3]; icol++ ){
								inp >> Z[ilay][irow][icol];
								cout << Z[ilay][irow][icol] << endl;
							}
						}
					}
				}
				else{
					datafile->getline(buffer,1024);
					//cout << buffer << endl;
					istringstream inp(buffer);
					for (int i = 0; i < Matrix_dim[3 - idim]; i++){
						inp >> X[idim][i];
						//cout << X[idim][i] << endl;
					}
				}
			}
		}

		{// read data for interpolation
			for (int ilay = 0; ilay < Matrix_dim[1]; ilay++){
				for (int irow = 0; irow < Matrix_dim[2]; irow++){
					datafile->getline(buffer,1024);
					istringstream inp(buffer);
					for (int icol = 0; icol < Matrix_dim[3]; icol++ ){
						inp >> V[ilay][irow][icol];
						//cout << V[ilay][irow][icol] << endl;
					}
				}
			}
		}
	}
}

template <int dim>
double interp<dim>::interpolate(double x0[]){
	int J[2]; int I[2]; int K[2];
	int nQ;
	double v0;
	double x[2];
	double y[2];
	double* z;
	double* Q;
	if (dim == 1){
		get_index(J, x, x0, 0);
		//cout << x[0] << " " << x[1] << endl;
		Q = new double[2]; nQ = 2;
		Q[0] = V[0][0][J[0]]; Q[1] = V[0][0][J[1]];
	}
	else if (dim == 2){
		get_index(J, x, x0, 0);
		//cout << x[0] << " " << x[1] << endl;
		get_index(I, y, x0, 1);
		//cout << y[0] << " " << y[1] << endl;
		Q = new double[4]; nQ = 4;
		Q[0] = V[0][I[0]][J[0]];	Q[1] = V[0][I[0]][J[1]];
		Q[2] = V[0][I[1]][J[0]];	Q[3] = V[0][I[1]][J[1]];
	}
	else if (dim == 3){
		get_index(J, x, x0, 0);
		//cout << x[0] << " " << x[1] << endl;
		get_index(I, y, x0, 1);
		//cout << y[0] << " " << y[1] << endl;
		if (LayElev == false){
			z = new double[2];
			get_index(K, z, x0, 2);
			//cout << z[0] << " " << z[1] << endl;
		}
		else{
			z = new double[10];
			get_index_Z(I, y, J, x, K, z, x0);
		}

		Q = new double[8]; nQ = 8;
		Q[0] = V[K[0]][I[0]][J[0]];	Q[1] = V[K[0]][I[0]][J[1]];
		Q[2] = V[K[0]][I[1]][J[0]];	Q[3] = V[K[0]][I[1]][J[1]];
		Q[4] = V[K[1]][I[0]][J[0]];	Q[5] = V[K[1]][I[0]][J[1]];
		Q[6] = V[K[1]][I[1]][J[0]]; Q[7] = V[K[1]][I[1]][J[1]];
	}
	//for (int i = 0; i < nQ; i++) cout << Q[i] << endl;

	int Check = CheckNodataValues(Q, nQ);

	if (Check == nQ){ // linear bilinear trilinear interpolation
		if (dim == 1){
			v0 = ( (Q[1] - Q[0]) / (x[1] - x[0]) ) * (x0[0] - x[0]) + Q[0];
		}
		else if (dim == 2){
			v0 = (  Q[0] * (x[1] - x0[0]) * (y[1] - x0[1]) +
					Q[1] * (x0[0] - x[0]) * (y[1] - x0[1]) +
					Q[2] * (x[1] - x0[0]) * (x0[1] - y[0]) +
					Q[3] * (x0[0] - x[0]) * (x0[1] - y[0]) ) / ( (x[1] - x[0]) * (y[1] - y[0]) );
		}
		else if (dim == 3 && LayElev == false){
			double xd, yd, zd;
			double c[4];
			double cc[2];
			xd = (x0[0] - x[0])/(x[1] - x[0]);
			yd = (x0[1] - y[0])/(y[1] - y[0]);
			zd = (x0[2] - z[0])/(z[1] - z[0]);

			c[0] = Q[0] * (1 - xd) + Q[1] * xd;
			c[1] = Q[2] * (1 - xd) + Q[3] * xd;
			c[2] = Q[4] * (1 - xd) + Q[5] * xd;
			c[3] = Q[6] * (1 - xd) + Q[7] * xd;

			cc[0] = c[0] * (1-yd) + c[1] * yd;
			cc[1] = c[2] * (1-yd) + c[3] * yd;

			v0 = cc[0] * (1 - zd) + cc[1] * zd;
		}
		else if (dim == 3 && LayElev == true){
			//interpolate the bottom and top layer using 2D interpolation
			// and then interpolate the actual value using 1D interpolation
			double topz, bottomz;
			interp<2> interp2D;
			double** tX;
			tX = new double*[2];
			tX[0] = new double [2];
			tX[1] = new double [2];
			double*** tV;
			tV = new double** [1];
			tV[0] = new double* [2];
			tV[0][0] = new double [2];
			tV[0][1] = new double [2];
			tX[0][0] = x[0]; tX[0][1] = x[1];
			tX[1][0] = y[0]; tX[1][1] = y[1];
			// bottom layer
			tV[0][0][0] = Q[0]; tV[0][0][1] = Q[1];
			tV[0][1][0] = Q[2]; tV[0][1][1] = Q[3];
			interp2D.get_data_c(2,2,1,tX,tV);
			bottomz = interp2D.interpolate(x0);

			// top layer
			tV[0][0][0] = Q[4]; tV[0][0][1] = Q[5];
			tV[0][1][0] = Q[6]; tV[0][1][1] = Q[7];
			interp2D.get_data_c(2,2,1,tX,tV);
			topz = interp2D.interpolate(x0);

			// 1D interpolation
			tX[0][0] = z[8]; tX[0][1] = z[9];
			tV[0][0][0] = bottomz; tV[0][0][1] = topz;
			interp<1> interp1D;
			interp1D.get_data_c(2,1,1,tX,tV);
			double z0[1];z0[0] = x0[2];
			v0 = interp1D.interpolate(z0);
		}
	}
	else{ //Inverse distance weighting
		double sum_w_u = 0;
		double sum_w = 0;
		double w;
		double **XX;
		XX = new double*[dim];
		if (Check <= 1){//find few nodes around x0 that have known values
			int **LIST;
			LIST = new int* [500];
			for (int i = 0; i < 500; i++){
				LIST[i] = new int[3];
				LIST[i][0] = -1; LIST[i][1] = -1; LIST[i][2] = -1;
			}
			if (dim == 1) {LIST[0][0] = J[0]; LIST[0][1] = 0; LIST[0][2] = 0;}
			else if (dim == 2 ){LIST[0][0] = J[0]; LIST[0][1] = I[0];LIST[0][2] = 0;}
			else if (dim == 3 ){LIST[0][0] = J[0]; LIST[0][1] = I[0]; LIST[0][2] = K[0];}

			int nlist = ListNeighbors(LIST);
			delete[] Q;
			Q = new double[nlist];
			nQ = nlist;
			for (int idim = 0; idim < dim; idim++)
				XX[idim] = new double[nlist];


			for (int i = 0; i < nlist; i++){
				Q[i] = V[ LIST[i][2] ][ LIST[i][1] ][ LIST[i][0] ];
				for (int idim = 0; idim < dim; idim++)
					XX[idim][i] = X[idim][LIST[i][idim]];
			}

			for (int i = 0; i < 500; i++)
				delete[] LIST[i];
			delete[] LIST;
		}
		else{
			if (dim == 2){
				XX[0] = new double[nQ];
				XX[1] = new double[nQ];
				XX[0][0] = x[0]; XX[0][1] = x[1]; XX[0][2] = x[0]; XX[0][3] = x[1];
				XX[1][0] = y[0]; XX[1][1] = y[1]; XX[1][2] = y[0]; XX[1][3] = y[1];
			}
			else if (dim == 3){
				XX[0] = new double[nQ];
				XX[1] = new double[nQ];
				XX[2] = new double[nQ];
				XX[0][0] = x[0]; XX[0][1] = x[1]; XX[0][2] = x[0]; XX[0][3] = x[1];
				XX[0][4] = x[0]; XX[0][5] = x[1]; XX[0][6] = x[0]; XX[0][7] = x[1];
				XX[1][1] = y[0]; XX[1][1] = y[0]; XX[1][2] = y[1]; XX[1][3] = y[1];
				XX[1][4] = y[0]; XX[1][5] = y[0]; XX[1][6] = y[1]; XX[1][7] = y[1];
				if (LayElev == false){
					XX[2][0] = z[0]; XX[2][1] = z[0]; XX[2][2] = z[0]; XX[2][3] = z[0];
					XX[2][4] = z[1]; XX[2][5] = z[1]; XX[2][6] = z[1]; XX[2][7] = z[1];
				}
				else if (LayElev == true){
					XX[2][0] = z[0]; XX[2][1] = z[1]; XX[2][2] = z[2]; XX[2][3] = z[3];
					XX[2][4] = z[4]; XX[2][5] = z[5]; XX[2][6] = z[6]; XX[2][7] = z[7];
				}
			}
		}

		double dst;
		for (int i = 0; i < nQ; i++){
			dst = 0;
			if (Q[i] == no_data)
				continue;
			for (int idim = 0; idim < dim; idim++){
				dst = dst + pow(x0[idim] - XX[idim][i],2);
			}
			dst = sqrt(dst);
			w = 1/pow(dst, 1.6);
			sum_w = sum_w + w;
			sum_w_u = sum_w_u + w * Q[i];
		}
		v0 = sum_w_u / sum_w;

		for (int i = 0; i < dim; i++)
			delete[] XX[i];
		delete[] XX;

	}

	delete[] Q;
	return v0;
}


template <int dim>
void interp<dim>::get_index(int I[], double x[], double x0[], int idim){
	int Nx = Matrix_dim[3 - idim];
	// index i means that x0 lays within the [i i+1] space
	// index -1 means that x0 < xmin
	//index Nx means that x0 > xmax
	if (x0[idim] < X[idim][0]){
		I[0] = 0; I[1] =0;
		x[0] = x0[idim] - (X[idim][1] - X[idim][0]);
		x[1] = X[idim][0];
	}
	else if (x0[idim] >= X[idim][Nx-1]){
		I[0] = Nx - 1; I[1] = I[0];
		x[0] = X[idim][Matrix_dim[3 - idim] - 1];
		x[1] = x0[idim] + (X[idim][Matrix_dim[3 - idim] - 1] - X[idim][Matrix_dim[3 - idim] - 2]);
	}
	else{
		for (int i = 0; i < Nx - 1; i++){
			if ((x0[idim] >= X[idim][i]) & (x0[idim] < X[idim][i+1])){
				I[0] = i;
				I[1] = i + 1;
				x[0] = X[idim][i];
				x[1] = X[idim][i + 1];
				break;
			}
		}
	}
}

template <int dim>
int interp<dim>::CheckNodataValues(double* Q, int nQ){
	int Ndata = 0;
	for (int i = 0; i < nQ; i++){
		if (Q[i] != no_data)
			Ndata++;
	}
	return Ndata;
}


template <int dim>
void interp<dim>::get_index_Z(int I[], double y[], int J[], double x[], int K[], double *z, double x0[]){
	// create a 2D interpolant to interpolate the z values
	double** tX;
	tX = new double*[2];
	tX[0] = new double [2];
	tX[1] = new double [2];
	double*** tV;
	tV = new double** [1];
	tV[0] = new double* [2];
	tV[0][0] = new double [2];
	tV[0][1] = new double [2];
	interp<2> tempinterp;
	double tempz;

	tX[0][0] = x[0]; tX[0][1] = x[1];
	tX[1][0] = y[0]; tX[1][1] = y[1];

	for (int k = 0; k <= Matrix_dim[1]; k++){
		if (k == 0){
			tV[0][0][0] = Z[k][I[0]][J[0]]; tV[0][0][1] = Z[k][I[0]][J[1]];
			tV[0][1][0] = Z[k][I[1]][J[0]]; tV[0][1][1] = Z[k][I[1]][J[1]];
			tempinterp.get_data_c(2,2,1,tX,tV);
			tempz = tempinterp.interpolate(x0);
			if (x0[2] < tempz){
				K[0] = 0; K[1] =0;
				z[0] = Z[0][I[0]][J[0]] - (Z[1][I[0]][J[0]] - Z[0][I[0]][J[0]]);
				z[1] = Z[0][I[0]][J[1]] - (Z[1][I[0]][J[1]] - Z[0][I[0]][J[1]]);
				z[2] = Z[0][I[1]][J[0]] - (Z[1][I[1]][J[0]] - Z[0][I[1]][J[0]]);
				z[3] = Z[0][I[1]][J[1]] - (Z[1][I[1]][J[1]] - Z[0][I[1]][J[1]]);
				z[4] = Z[0][I[0]][J[0]]; z[5] = Z[0][I[0]][J[1]];
				z[6] = Z[0][I[1]][J[0]]; z[7] = Z[0][I[1]][J[1]];
				z[8] = tempz - (Z[1][I[0]][J[0]] - Z[0][I[0]][J[0]]);
				z[9] = tempz;
				break;
			}
		}
		else if (k == Matrix_dim[1]){
			tV[0][0][0] = Z[k-1][I[0]][J[0]]; tV[0][0][1] = Z[k-1][I[0]][J[1]];
			tV[0][1][0] = Z[k-1][I[1]][J[0]]; tV[0][1][1] = Z[k-1][I[1]][J[1]];
			tempinterp.get_data_c(2,2,1,tX,tV);
			tempz = tempinterp.interpolate(x0);
			if (x0[2] > tempz){
				int Nz = Matrix_dim[1] - 1;
				K[0] = Nz; K[1] = Nz;
				z[0] = Z[Nz][I[0]][J[0]]; z[1] = Z[Nz][I[0]][J[1]];
				z[2] = Z[Nz][I[1]][J[0]]; z[3] = Z[Nz][I[1]][J[1]];
				z[4] = Z[Nz][I[0]][J[0]] + (Z[Nz][I[0]][J[0]] - Z[Nz-1][I[0]][J[0]]);
				z[5] = Z[Nz][I[0]][J[1]] + (Z[Nz][I[0]][J[1]] - Z[Nz-1][I[0]][J[1]]);
				z[6] = Z[Nz][I[1]][J[0]] + (Z[Nz][I[1]][J[0]] - Z[Nz-1][I[1]][J[0]]);
				z[7] = Z[Nz][I[1]][J[1]] + (Z[Nz][I[1]][J[1]] - Z[Nz-1][I[1]][J[1]]);
				z[8] = tempz;
				z[9] = tempz + (Z[Nz][I[0]][J[0]] - Z[Nz - 1][I[0]][J[0]]);
				break;
			}
		}
		else{
			tV[0][0][0] = Z[k][I[0]][J[0]]; tV[0][0][1] = Z[k][I[0]][J[1]];
			tV[0][1][0] = Z[k][I[1]][J[0]]; tV[0][1][1] = Z[k][I[1]][J[1]];
			tempinterp.get_data_c(2,2,1,tX,tV);
			double tempz1 = tempinterp.interpolate(x0);
			if (x0[2] >= tempz && x0[2] < tempz1){
				K[0] = k - 1; K[1] = k;
				z[0] = Z[k-1][I[0]][J[0]]; z[1] = Z[k-1][I[0]][J[1]];
				z[2] = Z[k-1][I[1]][J[0]]; z[3] = Z[k-1][I[1]][J[1]];
				z[4] = Z[k][I[0]][J[
				                    0]]; z[5] = Z[k][I[0]][J[1]];
				z[6] = Z[k][I[1]][J[0]]; z[7] = Z[k][I[1]][J[1]];
				z[8] = tempz;
				z[9] = tempz1;
				break;
			}
			tempz = tempz1;
		}
	}
}


template <int dim>
void interp<dim>::destroy(){
	for (int idim = 0; idim < dim; idim++){
		delete[] X[idim];
	}
	delete[] X;

	for (int ilay = 0; ilay < Matrix_dim[1]; ilay++){
		for (int irow = 0; irow < Matrix_dim[2]; irow++ ){
			delete[] V[ilay][irow];
			if (LayElev == true) delete[] Z[ilay][irow];
		}
		delete[] V[ilay];
		if (LayElev == true) delete[] Z[ilay];
	}
	delete[] V;
	if (LayElev == true) delete[] Z;

	cout << "Matrices destroyed!" << endl;

}


template <int dim>
void interp<dim>::set_LayElev(bool value){
	if (dim != 3)
		LayElev = false;
	else
		LayElev = value;
}

template <int dim>
int interp<dim>::ListNeighbors(int **&LIST){
	int istart = 0; int iend = 0;
	int ii, nFLT;
	int cnt_list = 1;
	int cnt_valid = 0;
	{	// check if the first value in the list is no_data
		int J, I, K;
		if (dim == 1) {J = LIST[0][0]; I = 0; K =0;}
		else if (dim == 2) {J = LIST[0][0]; I = LIST[0][1]; K = 0;}
		else if (dim == 3) {J = LIST[0][0]; I = LIST[0][1]; K = LIST[0][2];}
		if (V[K][I][J] != no_data)
		cnt_valid++;
	}
	int **FILTER;
	nFLT = Filter(FILTER, 1);
	//for (int iii = 0; iii < nFLT; iii++)
	//	cout << FILTER[iii][0] << " " << FILTER[iii][1] << endl;
	int Itest[3];
	int nGen = 10;
	bool exitloops = false;
	for (int i = 0; i < nGen; i++){
		for (int k = istart; k <= iend; k++){
			for (ii = 0; ii < nFLT; ii++){
				for (int idim = 0; idim < dim; idim++)
					Itest[idim] = LIST[k][idim] + FILTER[ii][idim];
				if (testneighbor(Itest, LIST)){
					int idim;
					for (idim = 0; idim < dim; idim++)
						LIST[cnt_list][idim] = Itest[idim];
					for (int ii = idim; ii < 3; ii++)
						LIST[cnt_list][ii] = 0;

					cnt_list++;
					{	// check if the new value in the list is no_data
						int J, I, K;
						if (dim == 1) {J = Itest[0]; I = 0; K =0;}
						else if (dim == 2) {J = Itest[0]; I = Itest[1]; K = 0;}
						else if (dim == 3) {J = Itest[0]; I = Itest[1]; K = Itest[2];}
						if (V[K][I][J] != no_data)
						cnt_valid++;
					}
					if (cnt_valid >= nP){
						exitloops = true;
						break;
					}
				}
				if (exitloops)
					break;
			}
			if (exitloops)
				break;
		}
		if (exitloops)
			break;
		istart = iend + 1;
		iend = cnt_list - 1;
	}

	nFLT = Filter(FILTER, 0);
	return cnt_list;
}

template <int dim>
bool interp<dim>::testneighbor(int *Itest, int **LIST){
	/**< Tests if the index is valid and exist in the #LIST.
	 * False if the index is not valid or exists in #LIST. Otherwise true
	*/
	bool outcome = true;
	for (int idim = 0; idim < dim; idim++){
		if (Itest[idim] < 0){
			outcome = false;
			return outcome;
		}
		if (Itest[idim] >= Matrix_dim[3 - idim]){
			outcome = false;
			return outcome;
		}
	}


	int J, I, K;
	if (dim == 1) {J = Itest[0]; I = 0; K =0;}
	else if (dim == 2) {J = Itest[0]; I = Itest[1]; K = 0;}
	else if (dim == 3) {J = Itest[0]; I = Itest[1]; K = Itest[2];}

	for (int j = 0; j < 500; j++){
		if (LIST[j][0] == -1)
			break;
		int tt = 0;
		for (int idim = 0; idim < dim; idim++)
			if (LIST[j][idim] == Itest[idim])
				tt++;
		if (tt == dim){
			outcome = false;
			return outcome;
		}
	}
	return outcome;
}

template <int dim>
int interp<dim>::Filter(int **&FLT, int mode){
	int cnt = 0;
	int nFLT;
	if (mode == 1){
		if (dim == 1){
			FLT = new int* [3];
			for (int j = 0; j < 3;j++)
				FLT[j] = new int[1];
			for (int j = -1; j < 2; j++){
				FLT[cnt][0] = j;
				cnt++;
			}
			nFLT = 3;
		}
		else if (dim == 2){
			FLT = new int* [9];
			for (int j = 0; j < 9;j++)
				FLT[j] = new int[2];
			for (int i = -1; i < 2; i++){
				for (int j = -1; j < 2; j++){
					FLT[cnt][0] = i;
					FLT[cnt][1] = j;
					cnt++;
				}
			}
			nFLT = 9;
		}
		else if (dim == 3){
			FLT = new int* [27];
			for (int j = 0; j < 27;j++)
				FLT[j] = new int[3];
			for (int k = -1; k < 2; k++){
				for (int i = -1; i < 2; i++){
					for (int j = -1; j < 2; j++){
						FLT[cnt][0] = i;
						FLT[cnt][1] = j;
						FLT[cnt][2] = k;
					}
				}
			}
			nFLT = 27;
		}
	}
	else if (mode == 0){
		if (dim == 1) nFLT = 3;
		else if (dim == 2) nFLT = 9;
		else if (dim == 3) nFLT = 27;

		for (int j = 0; j < nFLT;j++)
			delete[] FLT[j];
		delete[] FLT;
	}
	return nFLT;
}

template <int dim>
void interp<dim>::set_nP(int value){
	nP = value;
}

template <int dim>
void interp<dim>::set_Z(double ***tZ){
	if (dim == 3)
		Z = &tZ[0];
}
