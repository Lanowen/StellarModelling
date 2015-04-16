#include "Opal.hpp"
#include "Star.hpp"

Opal::Opal(Star* star, string filename, bool use_opal) : star(star), Graphable(2) {
	if (use_opal) {
		IFileHandle file(filename);

		table.resize(71);

		double temp;
		int row = 0;
		while (!file.eof()) {
			file >> temp;
			table[row].push_back(temp);
			if (file.peek() == '\n')
				row++;
		}
	}
}

double Opal::get_kappa(double temperature, double density) {
	if (!star->use_opal)
		return star->kappa.get();
	double T = log10(temperature);
	double R = log10(density / 1000.0 / pow(temperature * 1.0E-6, 3));

	int x = 1, y = 1;
	for (int i = 1; i < table.size(); i++) {
		if (table[i][0] <= T) {
			x = i;

		}
	}

	for (int i = 1; i < table[0].size(); i++) {
		if (table[0][i] <= R) {
			y = i;
		}
	}

	if (y == table[x].size() - 1 && (x == table.size() - 1 || y > table[x+1].size() - 1)) {
		return pow(10.0, table[x][y]) / 10.0; //divide by 10 to convert from cm�/g to m�/kg
	}
	if (x == table.size() - 1) {
		double fR = (R - table[0][y]) / (table[0][y + 1] - table[0][y]);
		double r1;

		r1 = table[x][y] + fR*(table[x][y + 1] - table[x][y]);

		return pow(10.0, r1) / 10.0; //divide by 10 to convert from cm�/g to m�/kg
	}
	if (y == table[x].size() - 1) {
		double fT = (T - table[x][0]) / (table[x + 1][0] - table[x][0]);
		double r1, r2;

		r1 = table[x][y];
		r2 = table[x + 1][y];

		return pow(10.0, (r1 + fT*(r2 - r1))) / 10.0; //divide by 10 to convert from cm�/g to m�/kg
	}


	double fT = (T - table[x][0]) / (table[x + 1][0] - table[x][0]), fR = (R - table[0][y]) / (table[0][y + 1] - table[0][y]);
	double r1, r2;
	
	r1 = table[x][y] + fR*(table[x][y + 1] - table[x][y]);
	r2 = table[x + 1][y] + fR*(table[x + 1][y + 1] - table[x + 1][y]);

	return pow(10.0, (r1 + fT*(r2 - r1)))/10.0; //divide by 10 to convert from cm�/g to m�/kg
}

void Opal::pushValues() {
	arr[0].push_back(log10(star->temperature.get()));
	arr[1].push_back(log10(get_kappa(star->temperature.get(),star->density.get())));
}