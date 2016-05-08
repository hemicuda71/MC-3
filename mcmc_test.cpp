#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string.h> // for memset
#include <time.h>
#include <stdlib.h>     /* srand, rand */
#include "timer.h"

using namespace std;

const double A1 = 4.15;
const double A2 = 0.531;
const double A3 = 67.3;

double cross_section(double E) {
	double core = A1 - A2*log(E);
	double tail = 1.0 - exp(-A3/E);
	return log(core*core * pow(tail,4.5));
}

// function to fit
// pos[0] = first coord, pos[1] = second coord, etc
// out = which output dim to return
// state = the parameters
double fit_func(double *pos, int out, double *state)
{
	// double Estar, E1 = state[2], g1 = state[0], g2 = state[1];
	// double E = pos[0];
	// Estar = pow(exp(cross_section(E1))/exp(cross_section(0.005)), 1.0/(g2-g1))
	//   * pow(E1, g2/(g2-g1)) / pow(0.005, g1/(g2-g1));
	// return (E < Estar ? cross_section(0.005) - g1*log(E/0.005)
	// 	              : cross_section(E1) - g2*log(E/E1));

	double Estar1, Estar2, Estar3, E = pos[0];
	double g1 = state[0], g2 = state[1], g3 = state[2], g4 = state[3];
	double E1 = state[4], E2 = state[5], E3 = state[6];
	Estar1 = pow(exp(cross_section(E1))/exp(cross_section(0.005)), 1.0/(g2-g1))
	  * pow(E1, g2/(g2-g1)) / pow(0.005, g1/(g2-g1));
	Estar2 = pow(exp(cross_section(E2))/exp(cross_section(E1)), 1.0/(g3-g2))
	  * pow(E2, g3/(g3-g2)) / pow(E1, g2/(g3-g2));
	Estar3 = pow(exp(cross_section(E3))/exp(cross_section(E2)), 1.0/(g4-g3))
	  * pow(E3, g4/(g4-g3)) / pow(E2, g3/(g4-g3));

	  return (E < Estar1 ? cross_section(0.005) - g1*log(E/0.005)
	  	   : (E < Estar2 ? cross_section(E1) - g2*log(E/E1)
	  	   : (E < Estar3 ? cross_section(E2) - g3*log(E/E2)
	  	   :               cross_section(E3) - g4*log(E/E3))));

	//return cross_section(0.005) - state[0] * log(pos[0]/0.005);
}

// fitness function (ie error function), typcially xi squared
// double *data = the data we are trying to fit. Organized for a given set of
//   coordinates of d_in dimensions, there will be d_out number of outputs
// nsize = number of positions
double fitness(double *pos, double *data, double *state, double uncertainty,
	int d_in, int d_out, int nsize)
{
	int i,j;
	double sum = 0.0, eval;

	// for all positions
	for(i = 0; i < nsize; i++) {
		// for each output
		for(j = 0; j < d_out; j++) {
			//cout << data[i*d_out + j] << ' ' << fit_func(pos + d_in*i, j, state) << endl;
			eval = data[i*d_out + j] - fit_func(pos + d_in*i, j, state);
			sum += eval*eval;
		}
	}
	//sum /= (uncertainty*uncertainty);

	return sum;
}

int main(void)
{
	srand(time(NULL));
	timer main_timer;
	int n, i, j, chain1, chain2, swap_count = 0;
	double STD_GUESS; // like the temperature of the search, 1 = full range, 0 = no range
	double EPS = 1.0, DELTA_TEMP, ANNEALING_FACTOR;
	double *xi_cur, *xi_prop, lhood_ratio, xi_min = 1.0E10, swap_percentage;

	// read set up params from file
	ifstream setupfile("setup.txt");
	int NSTATES, NCHAINS, NPARAMS, DIN, DOUT, NDATA, SWAP_RATE, PRINT_FREQ;
	setupfile >> NSTATES >> NCHAINS >> NPARAMS >> DIN >> DOUT >> NDATA >> DELTA_TEMP
	  >> ANNEALING_FACTOR >> STD_GUESS >> SWAP_RATE >> PRINT_FREQ;

	printf("NSTATES: %i  NCHAINS: %i NPARAMS: %i  DIN: %i  DOUT: %i  NDATA: %i\n",
		NSTATES, NCHAINS, NPARAMS, DIN, DOUT, NDATA);
	printf("DELTA_TEMP: %f  ANNEALING_FACTOR: %f  STD_GUESS: %f\n",
		DELTA_TEMP, ANNEALING_FACTOR, STD_GUESS);
	printf("SWAP_RATE: %i  PRINT_FREQ: %i\n", SWAP_RATE, PRINT_FREQ);

	double *state_proposed, *state_prev, *state_best, *heat, *state_t, *state_range;
	//state_array    = new double[NSTATES*NPARAMS];
	state_prev     = new double[NPARAMS*NCHAINS];
	state_proposed = new double[NPARAMS*NCHAINS];
	state_best     = new double[NPARAMS*NCHAINS];
	state_range    = new double[2*NPARAMS];
	state_t = new double[NPARAMS];
	heat    = new double[NCHAINS];
	xi_cur  = new double[NCHAINS];
	xi_prop = new double[NCHAINS];
	memset(state_prev, 0, NPARAMS*NCHAINS*sizeof(double));

	// set up heats for each chain
	// Note, we want DELTA_TEMP such that the swaps are accepted between 20-60% of the time
	for(j = 0; j < NCHAINS; j++) {
		heat[j] = 1.0 / (1.0 + DELTA_TEMP*j);
		printf("Heat %i: %f\n", j, heat[j]);
	}

	// set initial state
	for(i = 0; i < NPARAMS; i++) {
		setupfile >> state_prev[i] >> state_range[2*i] >> state_range[2*i+1];
		for(j = 1; j < NCHAINS; j++) // fill all other chains
			state_prev[i + j*NPARAMS] = (0.9 + 0.2*rand()/double(RAND_MAX))*state_prev[i];
		printf("Param %i: %f | range = [%f, %f]\n", i,
			state_prev[i], state_range[2*i], state_range[2*i+1]);
	}
	setupfile.close();

	// read position and data from file
	double *pos, *data;
	pos  = new double[NDATA*DIN];
	data = new double[NDATA*DOUT]; 
	ifstream posfile("pos.txt");
	ifstream datafile("data.txt");
	if(!posfile.good()) {
		cout << "pos.txt open fail, aborting program...\n\n";
		return -1;
	}
	if(!datafile.good()) {
		cout << "data.txt open fail, aborting program...\n\n";
		return -1;
	}
	
	for(i = 0; i < NDATA*DIN; i++){
		posfile >> pos[i];
		//cout << pos[i] << endl;
	}
	//cout << endl;
	for(i = 0; i < NDATA*DOUT; i++) {
		datafile >> data[i];
		//cout << data[i] << endl;
	}
	//cout << endl << endl;
	posfile.close();
	datafile.close();
	ofstream best_chain_file("best_chain.txt");
	ofstream chainfile("chain.txt");
	//////////////////////////////// MAINLOOP ////////////////////////////////////
	start_timer(main_timer);
	for(n = 0; n < NSTATES-1; n++) {
		// set proposed state for each chain
		for(i = 0; i < NCHAINS*NPARAMS; i++) {
			do {
				state_proposed[i] = state_prev[i]
				  + sqrt(state_range[2*(i%NPARAMS)+1] - state_range[2*(i%NPARAMS)])
				  * STD_GUESS*(2.0*rand()/double(RAND_MAX) - 1.0);
			} while((state_proposed[i] < state_range[2*(i%NPARAMS)] 
				||  state_proposed[i] > state_range[2*(i%NPARAMS)+1]));
		}
		// calc xi squared vals for current state and proposed state for each chain
		for(j = 0; j < NCHAINS; j++) {
			xi_cur[j]   = fitness(pos, data, state_prev     + j*NPARAMS, EPS, DIN, DOUT, NDATA);
			xi_prop[j]  = fitness(pos, data, state_proposed + j*NPARAMS, EPS, DIN, DOUT, NDATA);
			lhood_ratio = exp(heat[j]*(xi_cur[j] - xi_prop[j]));

			//cout << n << ' ' << lhood_ratio << ' ' << endl;
			if(lhood_ratio > rand()/double(RAND_MAX)) { // keep proposed state
				//cout << "new " << lhood_ratio << endl;
				if(j == 0 && (xi_prop[j] < xi_min)) { // only keep cold chain best
					xi_min = xi_prop[j];
					best_chain_file << n << ' ';
					for(i = 0; i < NPARAMS; i++) {
						best_chain_file << setprecision(9) << state_proposed[i] << ' ';
						state_best[i] = state_proposed[i];
					}
					best_chain_file << ' ' << xi_prop[j] << endl;
				}
				xi_cur[j] = xi_prop[j]; // save xi^2 to be used in swap calc later
				for(i = 0; i < NPARAMS; i++)
					state_prev[i + j*NPARAMS] = state_proposed[i + j*NPARAMS];
				
			} // else keep old state
		}
		// swap chains at random
		if(n%SWAP_RATE == 0) {
			chain1 = (NCHAINS-1.0) * rand()/double(RAND_MAX);
			do{
			chain2 = (NCHAINS-1.0) * rand()/double(RAND_MAX);
			} while(chain2 == chain1);
			lhood_ratio = exp((heat[chain1] - heat[chain2])*(xi_cur[chain1] - xi_cur[chain2]));
			//cout << "Swap liklihood: " << lhood_ratio << ' ' << chain1 << ' ' << chain2 << endl;
			if(lhood_ratio > rand()/double(RAND_MAX)) { // swap the chains
				swap_count++;
				//cout << swap_count << endl;
				for(i = 0; i < NPARAMS; i++) {
					state_t[i] = state_prev[i + chain1*NPARAMS];
					state_prev[i + chain1*NPARAMS] = state_prev[i + chain2*NPARAMS];
					state_prev[i + chain2*NPARAMS] = state_t[i];
				}
			}
		}

		// print current cold chain to file
		if(n%PRINT_FREQ == 0) {
			if(n%1000 == 0){
				STD_GUESS *= ANNEALING_FACTOR; // simulated annealing
				printf("Iter: %i  STD_GUESS = %f\n", n, STD_GUESS);
				current_time(main_timer);
				display_remaining_time(main_timer, n+1, NSTATES, 3);
			}
			chainfile << n << ' ' ;
			for(i = 0; i < NPARAMS; i++) // only print cold chain
				chainfile << setprecision(9) << state_prev[i] << ' ';
			// swap chain percentage
			swap_percentage = 100.0 * swap_count * SWAP_RATE / double(n+1.0);
			chainfile /*<< swap_percentage << ' '*/ << xi_cur[0] << endl;
		}
	}
	best_chain_file.close();
	chainfile.close();

	cout << "Best params: xi = " << xi_min << endl;
	for(i = 0; i < NPARAMS; i++) 
		cout << state_best[i] << ' ';
	cout << endl << endl;
	
	delete[] state_proposed, state_prev, state_best;
	delete[] pos, data, xi_cur, xi_prop, state_t, state_range;
	return 0;
}