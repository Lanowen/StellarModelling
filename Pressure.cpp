#include "Pressure.hpp"
#include "Star.hpp"
#include <iostream>

long double Pressure::get_P() {
	long double a, b, c;
	a = pow(3.0*pi*pi, 2.0 / 3.0)*pow(hbar, 2) / 5.0 / me*pow(star->density.get() / mp, 5.0 / 3.0);
	b = star->density.get()*kb*star->temperature.get() / star->mu.get_mu() / mp;
	c = 1.0 / 3.0*pow(star->temperature.get(), 4) * a_rad;
	return  (a+b+c);
}

long double Pressure::dP_T(long double rho, long double T) {
	return rho*kb / star->mu.get_mu() / mp + 4.0 / 3.0*pow(T, 3) * a_rad;
}

void Pressure::pushValues() {
	long double a, b, c;
	
	a = pow(3.0*pi*pi, 2.0 / 3.0)*pow(hbar,2) / 5.0 / me*pow(star->density.get() / mp, 5.0 / 3.0);
	b = star->density.get()*kb*star->temperature.get() / star->mu.get_mu() / mp;
	c = 1.0 / 3.0*pow(star->temperature.get(), 4) * a_rad;

	deg->push_back(a);
	therm->push_back(b);
	rad->push_back(c);
	arr[0].push_back(a+b+c);
}

long double Pressure::dP_rho(long double rho, long double T) {
	return pow(3.0*pi*pi, 2.0 / 3.0)*pow(hbar, 2) / 3.0 / me/mp*pow(rho / mp, 2.0 / 3.0) + kb*T / star->mu.get_mu() / mp;
}

long double Pressure::get_dP_T() {
	return dP_T(star->density.get(), star->temperature.get());
}

long double Pressure::get_dP_rho() {
	return dP_rho(star->density.get(), star->temperature.get());
}