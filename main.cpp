#include <iostream>
#include "Star.hpp"
#include <discpp.h>
#include <algorithm>
#include "IFileHandle.hpp"

using namespace std;



int main() {
	long double rho_c = 62900;
	Star* star = 0;
	double X = 0.7, Y = 0.28, Z = 0.02;
	long double t_c = 8.23E6L, rho_c_1 = rho_c, rho_c_2 = 0;
	bool ePP = true, eCNO = true, e3a = true;
	long double err_sensitivity = 1E-15L;

	long double shoot_delta_density = 200.0L;
	
	IFileHandle options_file("StarOptions.txt");
	ofstream stars_output("stars_tab_delimited.txt", ofstream::trunc);
	stars_output.precision(19);

	stars_output << "Star id" << "\t"
		<< "T_tau" << "\t"
		<< "L/Lsun" << "\t"
		<< "log10(T_tau)" << "\t"
		<< "M/Msun" << "\t"
		<< "R/Rsun" << "\t"
		<< "\t"
		<< "R/Rsun text" << "\t"
		<< "L/Lsun text" << "\t"
		<< "\t"
		<< "T_SB" << "\t"
		<< "log10(T_SB)" << endl;
	
	//for (int t = 16; t < 17; t++) {
		//long double t_c = 2.0E8, rho_c = 238407, rho_c_1 = rho_c, rho_c_2 = 0;
		//double X = 0.1, Y = 0.89, Z = 0.01;

		//long double t_c = 8.23E6L, rho_c_1 = rho_c, rho_c_2 = 0;
		
		//t_c = 5.0E6L + t*1E6L;
	string tag;
	bool graph_iteration = false;
	bool makepdf = true;
	bool use_opal = false;
	long double R_lim_default = 10.0L;

	while (!options_file.eof()) {
		
		int frac_test = 0;
		bool bisect = false;
		long double R_lim = R_lim_default;
		int t = 0;
		bool red_giant = false;
		long double He_cutoff = 0.0L;

		options_file >> tag;
		if (tag == "#" || tag.find("#") != tag.npos) {
			//do nothing
			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "shoot_density_change") {
			options_file >> shoot_delta_density;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "graph_every_iteration") {
			int temp;
			options_file >> temp;

			graph_iteration = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "make_pdf") {
			int temp;
			options_file >> temp;

			makepdf = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "XYZ") {
			options_file >> X >> Y >> Z;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "PP") {
			int temp;
			options_file >> temp;

			ePP = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "CNO") {
			int temp;
			options_file >> temp;

			eCNO = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "3a") {
			int temp;
			options_file >> temp;

			e3a = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "opal") {
			int temp;
			options_file >> temp;

			use_opal = temp;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "sensitivity") {
			options_file >> err_sensitivity;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "R_lim_default") {
			options_file >> R_lim_default;

			R_lim = R_lim_default;

			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
		}
		else if (tag == "star") {
			options_file >> t >> rho_c >> t_c;
			red_giant = false;
		}
		else if (tag == "redgiant") {
			options_file >> t >> rho_c >> t_c >> He_cutoff;
			red_giant = true;
		}
		else {
			//do nothing
			options_file.ignore(numeric_limits<streamsize>::max(), '\n');
			continue;
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
			star = new Star(t_c, rho_c, X, Y, Z, 1.0L, LDBL_MAX, 1.2L*R_lim + 1.0L, ePP, eCNO, e3a, use_opal, err_sensitivity, "opal_.7_.28_.02.txt", He_cutoff);
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
			int count_min = 60;
			while (true){
				rho_c = (rho_c_1 + rho_c_2) / 2.0;
				if (star != 0)
					delete star;
				star = new Star(t_c, rho_c, X, Y, Z, 1.0L, LDBL_MAX, 1.2L*R_lim + 1.0L, ePP, eCNO, e3a, use_opal, err_sensitivity, "opal_.7_.28_.02.txt", He_cutoff);

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

				long double frac = star->frac_diff();

				if (last_pos_frac == 0 || abs(last_pos_frac->frac_diff()) > abs(frac)) {
					if (last_pos_frac != 0)
						delete last_pos_frac;
					last_pos_frac = star;
					star = 0;
				}

				if (frac > 0) { //positive, go lower			 
				down:
					lastDir = 2;
					R_lim = R_lim;
					cout << endl << endl << "Bisecting down." << endl;
					rho_c_1 = min(rho_c_1, rho_c_2);
					rho_c_2 = rho_c;
				}
				else{ //negative, go higher
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
		stars_output << t << "\t"
			<< last_pos_frac->T_star << "\t"
			<< last_pos_frac->L_max / Lsun << "\t"
			<< log10(last_pos_frac->T_star) << "\t"
			<< last_pos_frac->M_max / Msun << "\t"
			<< last_pos_frac->R_star / Rsun << "\t"
			<< "\t"
			<< ((last_pos_frac->M_max / Msun < 1.66) ? (1.06*pow(last_pos_frac->M_max / Msun, 0.945)) : (1.33*pow(last_pos_frac->M_max / Msun, 0.555))) << "\t"
			<< ((last_pos_frac->M_max / Msun < 0.7) ? (0.35*pow(last_pos_frac->M_max / Msun, 2.62)) : (1.02*pow(last_pos_frac->M_max / Msun, 3.92))) << "\t"
			<< "\t"
			<< pow(last_pos_frac->L_max / 4.0 / pi / pow(last_pos_frac->R_star, 2) / sigma_B, 1.0 / 4.0) << "\t"
			<< log10(pow(last_pos_frac->L_max / 4.0 / pi / pow(last_pos_frac->R_star, 2) / sigma_B, 1.0 / 4.0)) << endl;

		delete last_pos_frac;

		options_file.ignore(numeric_limits<streamsize>::max(), '\n');
	}

	stars_output.close();

	system("PAUSE");

	return 0;
}