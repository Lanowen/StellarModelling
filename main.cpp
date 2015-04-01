#include <iostream>
#include "Star.hpp"
#include <discpp.h>
#include <algorithm>
#include "IFileHandle.hpp"

using namespace std;

int main() {
	//233840764
	//long double rho_c = 58560;
	long double rho_c = 62900;
	//long double rho_c = 83480;
	Star* star = 0;
	double X = 0.7, Y = 0.28, Z = 0.02;
	long double t_c = 8.23E6L, rho_c_1 = rho_c, rho_c_2 = 0;
	bool ePP = true, eCNO = true, e3a = true;
	long double err_sensitivity = 1E-15L;

	long double shoot_delta_density = 200.0L;
	
	
	IFileHandle opt("StarOptions.txt");
	
	//for (int t = 16; t < 17; t++) {
		//long double t_c = 2.0E8, rho_c = 238407, rho_c_1 = rho_c, rho_c_2 = 0;
		//double X = 0.1, Y = 0.89, Z = 0.01;

		//long double t_c = 8.23E6L, rho_c_1 = rho_c, rho_c_2 = 0;
		
		//t_c = 5.0E6L + t*1E6L;
	string tag;
	bool graph_iteration = false;
	bool makepdf = true;
	long double R_lim_default = 10.0L;

	while (!opt.eof()) {
		
		int frac_test = 0;
		bool bisect = false;
		long double R_lim = R_lim_default;
		int t = 0;
		bool red_giant = false;
		long double He_cutoff = 0.0L;

		opt >> tag;
		if (tag == "#" || tag.find("#") != tag.npos) {
			//do nothing
			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "shoot_density_change") {
			opt >> shoot_delta_density;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "graph_every_iteration") {
			int temp;
			opt >> temp;

			graph_iteration = temp;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "make_pdf") {
			int temp;
			opt >> temp;

			makepdf = temp;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "XYZ") {
			opt >> X >> Y >> Z;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "PP") {
			int temp;
			opt >> temp;

			ePP = temp;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "CNO") {
			int temp;
			opt >> temp;

			eCNO = temp;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "3a") {
			int temp;
			opt >> temp;

			e3a = temp;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "sensitivity") {
			opt >> err_sensitivity;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "R_lim_default") {
			opt >> R_lim_default;

			R_lim = R_lim_default;

			opt.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "star") {
			opt >> t >> rho_c >> t_c;
			red_giant = false;
		}
		else if (tag == "redgiant") {
			opt >> t >> rho_c >> t_c >> He_cutoff;
			red_giant = true;
		}

		//t--;
		rho_c_1 = rho_c;
		rho_c_2 = 0;				

		std::streamsize ss = std::cout.precision();		

		for (int i = 0; i < 10000; i++) {
			if (star != 0)
				delete star;
			std::cout.precision(19);
			cout << endl << endl << i << " /10000" << "=== Testing density: " << rho_c << " kg/m^3 , T = " << t_c << " K ===" << endl << endl;
			std::cout.precision(ss);
			star = new Star(t_c, rho_c, X, Y, Z, 1.0L, LDBL_MAX, 1.2L*R_lim + 1.0L, ePP, eCNO, e3a, err_sensitivity, He_cutoff);
			RK4::step = 5000;
			star->solve();

			if (graph_iteration) {
				star->graph(t, makepdf);
			}

			double frac = star->frac_diff();

			if (frac > 0 && !isinf(frac)) { //positive, go lower
				R_lim = star->R_star / Rsun;
				if (frac_test == 1 && !isnan(frac) && !isinf(frac)) {
					bisect = true;
					cout << endl << endl << "===Beginning bisection===" << endl;
					break;
				}
				else {
					if (!isnan(frac) && !isinf(frac))
						frac_test = 2;
					rho_c_1 = rho_c;
					if (rho_c - shoot_delta_density < 0)
						rho_c /= 2;
					else 
						rho_c -= shoot_delta_density;
				}
			}
			else{ //negative, go higher
				if (frac_test == 2 && !isnan(frac) && !isinf(frac)) {
					bisect = true;
					cout << endl << endl << "===Beginning bisection===" << endl;
					break;
				}
				else {
					if (!isnan(frac) && !isinf(frac))
						frac_test = 1;
					rho_c_1 = rho_c;
					rho_c += shoot_delta_density;
				}
			}
		}


		Star* last_pos_frac = 0;

		if (bisect) {
			rho_c_2 = rho_c;
			int lastDir = 0;
			int count = 0;
			int count_min = 65;
			while (true){
				rho_c = (rho_c_1 + rho_c_2) / 2.0;
				if (star != 0)
					delete star;
				star = new Star(t_c, rho_c, X, Y, Z, 1.0L, LDBL_MAX, 1.2L*R_lim + 1.0L, ePP, eCNO, e3a, err_sensitivity, He_cutoff);

				std::cout.precision(19);
				cout << count << " /" << count_min << (count > count_min ? " (trying to end on positive fractional) " : "")  << "=== Testing density: " << rho_c << " kg/m^3 , T = " << t_c << " K ===" << endl << endl;
				std::cout.precision(ss);

				star->solve();

				if (graph_iteration) {
					star->graph(t, makepdf);
				}

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

				if (star->frac_diff() > 0) {
					if (last_pos_frac == 0 || last_pos_frac->frac_diff() > star->frac_diff()) {
						if (last_pos_frac != 0)
							delete last_pos_frac;
						last_pos_frac = star;
						star = 0;
					}					
				down:
					lastDir = 2;
					R_lim = R_lim;
					cout << endl << endl << "Bisecting down." << endl;
					rho_c_1 = min(rho_c_1, rho_c_2);
					rho_c_2 = rho_c;
				}
				else{
				up:

					lastDir = 1;
					cout << endl << endl << "Bisecting up." << endl;
					rho_c_1 = max(rho_c_1, rho_c_2);
					rho_c_2 = rho_c;
				}

				if (count >= count_min)
					break;

				count++;
			}
		}

		last_pos_frac->graph(t, makepdf);
		delete last_pos_frac;

		opt.ignore(numeric_limits<streamsize>::max(), '\n');
	}

	//Star* star;
	//star = new Star(8.23E6, 58560, 0.70, 0.28, 0.02, 5000.0L, 10000.0L, 1.2L);
	//star = new Star(8.23E6, 58546.296424781904, 0.70, 0.28, 0.02);
	//star = new Star(8.23E6, 57001.40863806593, 0.70, 0.28, 0.02);
	//star = new Star(1.571E7, 1.622E5, 0.70, 0.28, 0.02);

	//star->solve();
	//star->graph(1);

	system("PAUSE");

	return 0;
}