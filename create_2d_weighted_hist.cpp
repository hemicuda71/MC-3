// 2d hist with weighted values

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string.h> // for memset
#include <stdlib.h>     /* atoi */

using namespace std;

const int XRES = 100;
const int YRES = 100;
const double BURN_IN = 0.3;

int main(int argc, char* argv[])
{
	int i, j, n, XPARAM, YPARAM;
	ifstream setupfile("setup.txt");
	ifstream chain_file("chain.txt");

	// parameters to plot together
	XPARAM = atoi(argv[1]);
	YPARAM = atoi(argv[2]);

	// read setup file
	double DELTA_TEMP, ANNEALING_FACTOR, STD_GUESS;
	int NSTATES, NCHAINS, NPARAMS, DIN, DOUT, NDATA, SWAP_RATE, PRINT_FREQ;
	setupfile >> NSTATES >> NCHAINS >> NPARAMS >> DIN >> DOUT >> NDATA >> DELTA_TEMP
	  >> ANNEALING_FACTOR >> STD_GUESS >> SWAP_RATE >> PRINT_FREQ;

	printf("NSTATES: %i  NCHAINS: %i NPARAMS: %i  DIN: %i  DOUT: %i  NDATA: %i\n",
		NSTATES, NCHAINS, NPARAMS, DIN, DOUT, NDATA);
	printf("DELTA_TEMP: %f  ANNEALING_FACTOR: %f  STD_GUESS: %f\n",
		DELTA_TEMP, ANNEALING_FACTOR, STD_GUESS);
	printf("SWAP_RATE: %i  PRINT_FREQ: %i\n", SWAP_RATE, PRINT_FREQ);
	cout << "BURN_IN fraction: " << BURN_IN << endl;

	string param_names[] = {"powerlaw1", "powerlaw2", "powerlaw3", "powerlaw4",
		"E1", "E2", "E3", ".txt"};
	string h_string("2dhist_");
	string file_name(h_string + param_names[XPARAM] + param_names[YPARAM] + param_names[NPARAMS]);
	char file_name_char[file_name.size()+1];
	strcpy(file_name_char, file_name.c_str());
	//const char *param_names[] = {"powerlaw1", "powerlaw2", "powerlaw3", "powerlaw4",
	//	"E1", "E2", "E3"};
	double *state_range, *state_cur, dud, xi;
	state_range = new double[2*NPARAMS];
	state_cur   = new double[NPARAMS];

	for(i = 0; i < NPARAMS; i++)
		setupfile >> dud >> state_range[2*i] >> state_range[2*i+1];
	setupfile.close();

	double XMIN = state_range[2*XPARAM], XMAX = state_range[2*XPARAM + 1];
	double YMIN = state_range[2*YPARAM], YMAX = state_range[2*YPARAM + 1];

	int xid, yid;
	double *histogram;
	const int HIST_SIZE = XRES * YRES;
	histogram = new double[HIST_SIZE];
	memset(histogram, 0, HIST_SIZE*sizeof(double));

	// read chain file
	for(n = 0; n < NSTATES/PRINT_FREQ; n++) {
		chain_file >> dud;
		for(i = 0; i < NPARAMS; i++)
			chain_file >> state_cur[i];
		chain_file >> xi;
		if(n > BURN_IN * NSTATES/PRINT_FREQ) {
			xid = (XRES-1) * (state_cur[XPARAM] - XMIN) / (XMAX - XMIN);
			yid = (YRES-1) * (state_cur[YPARAM] - YMIN) / (YMAX - YMIN);
			histogram[xid + yid*XRES] += 1.0 / xi;
		}
	}

	// print to file to be plotted in SciDavis
	ofstream hist_file(file_name_char);
	cout << "Printing histogram file: " << file_name_char << endl << endl;
	for(j = 0; j < YRES; j++) {
		for(i = 0; i < XRES; i++) {
			hist_file << (histogram[i + j*XRES] > 0 ? log(histogram[i + j*XRES]) : 0.0) << ' ' ;
		}
		hist_file << endl;
	}

	hist_file.close();

	delete[] state_range, state_cur, histogram;
	return 0;
}