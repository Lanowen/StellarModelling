#include <iostream>
#include "Star.hpp"
#include <discpp.h>
#include <algorithm>

using namespace std;

int main() {
	long double t_c = 8.23E6, rho_c = 58560, rho_c_1 = rho_c, rho_c_2 = 0;
	//long double t_c = 5.23E6, rho_c = 60000, rho_c_1 = rho_c, rho_c_2 = 0;
	
	Star* star = 0;
	RK4::step = 10000;
	int frac_test = 0;
	bool bisect = false;

	std::streamsize ss = std::cout.precision();

	for (int i = 0; i < 100; i++) {
		if (star == 0)
			delete star;
		std::cout.precision(19);
		cout << endl << endl << i << " /100" << "=== Testing density: " << rho_c << " kg/m^3 ===" << endl << endl;
		std::cout.precision(ss);
		star = new Star(t_c, rho_c, 0.70, 0.28, 0.02);

		star->solve();

		double frac = star->frac_diff();
		
		if (frac < 0) { //negative go higher
			if (frac_test == 2 && !isnan(frac) && !isinf(frac)) {
				bisect = true;
				cout << endl << endl << "===Beginning bisection===" << endl;
				break;
			}
			else {
				if (!isnan(frac) && !isinf(frac))
					frac_test = 1;
				rho_c_1 = rho_c;
				rho_c += 100;
			}
		}
		else { //positive go lower
			if (frac_test == 1 && !isnan(frac) && !isinf(frac)) {
				bisect = true;
				cout << endl << endl << "===Beginning bisection===" << endl;
				break;
			}
			else {
				if (!isnan(frac) && !isinf(frac))
					frac_test = 2;
				rho_c_1 = rho_c;
				rho_c -= 100;
			}
		}	
	}

	if (bisect) {
		rho_c_2 = rho_c;
		int lastDir = 0;
		for (int i = 0; i < 50; i++) {
			rho_c = (rho_c_1 + rho_c_2) / 2.0;
			delete star;
			star = new Star(t_c, rho_c, 0.70, 0.28, 0.02);

			std::cout.precision(19);
			cout << i << " /50" << "=== Testing density: " << rho_c << " kg/m^3 ===" << endl << endl;
			std::cout.precision(ss);

			star->solve();

			if (isnan(star->frac_diff()) || isinf(star->frac_diff())) {
				if (lastDir == 1)
					goto down;
				else if (lastDir == 2)
					goto up;
				else {
					cout << "Some error" << endl;
					break;
				}
			}

			if (star->frac_diff() < 0) {
				up:
				lastDir = 1;
				cout << endl << endl << "Bisecting up." << endl;
				rho_c_1 = max(rho_c_1, rho_c_2);
				rho_c_2 = rho_c;
			}
			else {
				down:
				lastDir = 2;
				cout << endl << endl << "Bisecting down." << endl;
				rho_c_1 = min(rho_c_1, rho_c_2);
				rho_c_2 = rho_c;
			}
		}
	}

	star->graph();

	//star = new Star(8.23E6, 58560, 0.70, 0.28, 0.02);
	//star = new Star(8.23E6, 58546.296424781904, 0.70, 0.28, 0.02);
	//star = new Star(1.571E7, 1.622E5, 0.70, 0.28, 0.02);

	//star->solve();
	//star->graph();

	system("PAUSE");

	return 0;
}