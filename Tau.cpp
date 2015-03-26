#include "Tau.hpp"
#include "Star.hpp"

long double Tau::dtau_r(long double r, long double tau){
	//return star->opal.get_kappa(star->temperature.get(), star->density.get()) * star->density.get();
	return star->kappa.get() * star->density.get();
}