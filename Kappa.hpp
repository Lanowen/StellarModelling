#pragma once

#include "Graphable.hpp"
#include "Mu.hpp"
#include <algorithm>
#include "Star_Constants.hpp"
#include <vector>

class Star;

class Kappa : public Graphable {
public:
	vector<double> *kes, *kff, *kh, *k, *logT;
	Star* star;

	long double get();

	Kappa(Star* star) : star(star), Graphable(5) {
		kes = &arr[0];
		kff = &arr[1];
		kh = &arr[2];
		k = &arr[3];
		logT = &arr[4];
	}

	virtual void pushValues();
};