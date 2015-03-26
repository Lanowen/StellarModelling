#include "Kappa.hpp"
#include "Star.hpp"
#include <iostream>

#define _kh (2.5E-32*(star->mu.get_Z() / 0.02)*pow(star->density.get()/1E3, 0.5)*pow(star->temperature.get(), 9))
#define _kes (0.02*(1.0 + star->mu.get_X()))
#define _kff (1.0E24*(1.0 + star->mu.get_X())*(star->mu.get_Z() + 0.0001)*pow(star->density.get() / 1E3, 0.7)*pow(star->temperature.get(), -3.5))

long double Kappa::get() {
	return 1.0 / (1.0 / _kh + 1.0 / max(_kes, _kff));
}

void Kappa::pushValues() {
	kh->push_back(_kh);
	kes->push_back(_kes);
	kff->push_back(_kff);
	logT->push_back(log10(star->temperature.get()));
	k->push_back(1.0 / (1.0 / _kh + 1.0 / max(_kes, _kff)));
}