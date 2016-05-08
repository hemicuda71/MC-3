#include <iostream>
#include <fstream>
#include <string.h> // for memset

using namespace std;

// Maximum storage to allocate
const size_t CHUNK_SIZE = 3E7;

// Resolution of 3d Cartesian grid
const int X_RES = 100;
const int Y_RES = 100;
const int Z_RES = 100;

// Physical dimensions of bounding box of data
const float X_MIN = 0.02;
const float X_MAX = 0.3;
const float Y_MIN = 3.0;
const float Y_MAX = 5.0;
const float Z_MIN = 0.5;
const float Z_MAX = 15.0;

// i = x, j = y, k = z
// nn = [nx, ny, nz]
int arrayIdx3(int i, int j, int k, int* nn) 
{
	return i + nn[0]*(j + nn[1]*k);
}

int main(void) {
	
	const char data_file_name[] = "chain.txt";
	const char vtk_file_name[]  = "chain_1_4_5.vtk";
	int count = 0, i, j, k, x_id, y_id, z_id, tot_chunks = 0;
	int nn[3] = {X_RES, Y_RES, Z_RES};

	// set up histogram
	const size_t HIST_SIZE = X_RES * Y_RES * Z_RES;
	cout << "HIST_SIZE = " << HIST_SIZE << endl;
	float *histogram;
	histogram = new float[HIST_SIZE];
	memset(histogram, 0, HIST_SIZE*sizeof(float));	

	float *data_chunk, dud, *weight_chunk;
	data_chunk   = new float[3*CHUNK_SIZE];
	weight_chunk = new float[CHUNK_SIZE];

	// open data file
	ifstream data_file(data_file_name);
	cout << endl << data_file_name << endl;
	if(data_file.good()) {
		cout << " File open successfully\n";
	} else {
		cout << " File open failed\n";
		return -1;
	}

	/// Algorithm
	// 1) Read data until RAM_SIZE amount of data read into array
	// 2) Sort the data
	// 3) Bin the chunk of data with resolution X_RES*Y_RES*Z_RES
	// 4) Repeat until file is finished
	// 5) Print 3d histogram data to vtk file
	while(!data_file.eof()) {
		// 1) fill chunk with data
		count = 0;
		while(!data_file.eof() && count < CHUNK_SIZE) {
			// data_file >> data_chunk[3*count]
			// 	  	  >> data_chunk[1 + 3*count]
			// 	  	  >> data_chunk[2 + 3*count];
			data_file >> dud >> data_chunk[3*count] >> dud >> dud
					  >> data_chunk[3*count+1] >>  data_chunk[3*count+2] >> dud
					  >> dud >> weight_chunk[count];
			count++;
		}
		tot_chunks++;
		printf("Chunk %i read complete\n", tot_chunks);

		// 2) sort data

		// 3) Bin data
		for(i = 0; i < count; i++) {
			x_id = (X_RES-1) * (data_chunk[3*i]     - X_MIN) / (X_MAX - X_MIN);
			y_id = (Y_RES-1) * (data_chunk[3*i + 1] - Y_MIN) / (Y_MAX - Y_MIN);
			z_id = (Z_RES-1) * (data_chunk[3*i + 2] - Z_MIN) / (Z_MAX - Z_MIN);

			// cout << arrayIdx3(x_id, y_id, z_id, nn) << '/' << HIST_SIZE;
			// cout << '\t' << x_id << ' ' << y_id << ' ' << z_id << endl;

			histogram[arrayIdx3(x_id, y_id, z_id, nn)] += 1.0/weight_chunk[i];
		}
		
	}
	data_file.close();
	cout << "Finished binning all data\n";

	// 5) Print vtk file
	float x_pos, y_pos, z_pos;
	int POINTS = (X_RES+1)*(Y_RES+1)*(Z_RES+1);
	int CELLS  =  HIST_SIZE;

	ofstream vtk_file(vtk_file_name);
	vtk_file << "# vtk DataFile Version 3.0\n";
	vtk_file << "3d point to 3d histrogram\n";
	vtk_file << "ASCII\n";
	vtk_file << "DATASET RECTILINEAR_GRID\n";
	vtk_file << "DIMENSIONS " << X_RES+1 << ' ' << Y_RES+1 << ' ' << Z_RES+1 << endl;
	vtk_file << "X_COORDINATES " << X_RES+1 << " float\n";
	for(i = 0; i <= X_RES; i++)
		vtk_file << X_MIN + (X_MAX - X_MIN) * i / float(X_RES) << endl;
	vtk_file << "Y_COORDINATES " << Y_RES+1 << " float\n";
	for(j = 0; j <= Y_RES; j++)
		vtk_file <<  Y_MIN + (Y_MAX - Y_MIN) * j / float(Y_RES) << endl;
	vtk_file << "Z_COORDINATES " << Z_RES+1 << " float\n";
	for(k = 0; k <= Z_RES; k++)
		vtk_file << Z_MIN + (Z_MAX - Z_MIN) * k / float(Z_RES) << endl;

	vtk_file << "\nCELL_DATA " << CELLS << endl;
	vtk_file << "FIELD FieldData 1\n";
	vtk_file << "bin_data 1 " << CELLS << " float\n";
	for(i = 0; i < CELLS; i++)
		vtk_file << histogram[i] << endl;

	cout << "Finished printing VTK file\n\n";

	delete[] histogram, data_chunk, weight_chunk;
	return 0;
}