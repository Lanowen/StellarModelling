#include "Opal.hpp"
#include "Star.hpp"

double Opal::get_kappa(double temperature, double density) {
	//return star->kappa.get();
	double T = log10(temperature);
	double R = log10(density / 1000.0 / pow(temperature / 1E6, 3));

	int x = 0, y = 0;
	for (int i = 0; i < table.size(); i++) {
		if (table[i][0] < T) {
			x = i;

		}
	}

	for (int i = 0; i < table[0].size(); i++) {
		if (table[0][i] < R) {
			y = i;
		}
	}

	if (x == table.size()-1) {
		stringstream ss;
		ss << "Out of bounds of OPAL table values, cannot extrapolate data. temperature = " << temperature << ". logT = " << T << ". density = " << density << ". logR = " << R << ".";
		throw OpalException(ss.str().c_str());
		return pow(10.0, 9.9);
		//return star->kappa.get();
	}
	if (y == table[x].size()-1) {
		stringstream ss;
		ss << "Out of bounds of OPAL table values, cannot extrapolate data. density = " << density << ". logR = " << R << ". temperature = " << temperature << ". logT = " << T << ".";
		throw OpalException(ss.str().c_str());
		return pow(10.0, 9.9);
		//return star->kappa.get();
	}

	double fT = (T - table[x][0]) / (table[x + 1][0] - table[x][0]), fR = (R - table[0][y]) / (table[0][y + 1] - table[0][y]);
	double r1, r2;
	r1 = table[x][y] + fR*(table[x][y + 1] - table[x][y]);
	r2 = table[x + 1][y] + fR*(table[x + 1][y + 1] - table[x + 1][y]);

	return pow(10.0, (r1 + fT*(r2 - r1)))/10.0; //divide by 10 to convert from cm²/g to m²/kg
}

void Opal::pushValues() {
	arr[0].push_back(log10(star->temperature.get()));
	arr[1].push_back(get_kappa(star->temperature.get(),star->density.get()));
}