#include "Energy.hpp"
#include "Star.hpp"

long double Energy::get() {
	long double e_pp = 1.07E-7 * star->density.get() / 1E5 * pow(star->mu.get_X(), 2) * pow(star->temperature.get() / 1E6, 4);
	long double e_CNO = 8.24E-26 * star->density.get() / 1E5 * 0.03 * pow(star->mu.get_X(), 2) * pow(star->temperature.get() / 1E6, 19.9);
	long double e_3alpha = 3.85E-8 * pow(star->density.get() / 1E5, 2) * pow(star->mu.get_Y(), 3) * pow(star->temperature.get() / 1E8, 44);
	return e_pp + e_CNO;// +e_3alpha;
}