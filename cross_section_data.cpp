// g++ -std=c++11 -O3 cross_section_data.cpp

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

const double A1 = 4.15;
const double A2 = 0.531;
const double A3 = 67.3;

double cross_section(double E) {
	double core = A1 - A2*log(E);
	double tail = 1.0 - exp(-A3/E);
	return log(core*core * pow(tail,4.5));
}

int main(int argc, char* argv[])
{
	double Emin = 0.005, Emax = 250.0, E;
	string str(argv[1]);
	size_t sz;
	int Ne = std::stoi(str, &sz); //= argv[1];
	//cout << Ne << endl;

	ofstream datafile("data.txt");
	ofstream posfile("pos.txt");

	for(int i = 0; i < Ne; i++) {
		E = Emin + (Emax-Emin) * i / double(Ne-1.0);
		posfile << setprecision(10) << E << endl;
		datafile << setprecision(10) << cross_section(E) << endl;
	}
	datafile.close();
	posfile.close();

	return 0;
}
