#include "Temperature.hpp"
#include "Star.hpp"
#include <iostream>

long double Temperature::dT_dr_conv(long double r, long double T) {
	return (1.0 - (1.0 / gamma))*(T / star->pressure.get_P())*G*star->mass.get()*star->density.get() / pow(r, 2) ;
}

long double Temperature::dT_dr_rad(long double r, long double T) {
	//return 3.0*star->opal.get_kappa(T, star->density.get())*star->density.get()*star->luminosity.get() / (pow(T, 3) * 16.0 * a_rad * c_0 * pi * pow(r, 2));
	return 3.0*star->kappa.get()*star->density.get()*star->luminosity.get() / (pow(T, 3) * 16.0 * a_rad * c_0 * pi * pow(r, 2));
}

long double Temperature::dT_dr(long double r, long double T) {
	long double conv = dT_dr_conv(r, T);
	long double rad = dT_dr_rad(r, T);
	return -min(dT_dr_conv(r, T), dT_dr_rad(r, T)) ;
}

long double Temperature::criterion() { //dlnTdlnP
	return (3.0*star->kappa.get()*star->luminosity.get()*star->pressure.get_P()) / (pow(get(), 4) * 16.0 * a_rad * c_0 * pi *G * star->mass.get());
	//return (3.0*star->opal.get_kappa(get(), star->density.get())*star->luminosity.get()*star->pressure.get_P()) / (pow(get(), 4) * 16.0 * a_rad * c_0 * pi *G * star->mass.get());
}