#pragma once

#include <algorithm>
#include "Graphable.hpp"

class Star;

#define _e_pp (1.07E-7 * star->density.get() / 1E5 * pow(star->mu.get_X(), 2) * pow(star->temperature.get() / 1E6, 4))
#define _e_CNO (8.24E-26 * star->density.get() / 1E5 * 0.03 * pow(star->mu.get_X(), 2) * pow(star->temperature.get() / 1E6, 19.9))
#define _e_3alpha (3.85E-8 * pow(star->density.get() / 1E5, 2) * pow(star->mu.get_Y(), 3) * pow(star->temperature.get() / 1E8, 44))
//#define _e_3alpha (2.9E5 * pow(star->density.get() / 1E3, 2) * pow(star->mu.get_Y(), 3)/pow(star->temperature.get()/1E9, 3)*exp(-4.4109/(star->temperature.get()/1E9)))

class Energy : public Graphable{
public:
	Star* star;

	long double get();

	Energy(Star* star) : star(star), Graphable(4) {}

	virtual void pushValues();

	inline virtual void popValues() {
		arr[0].pop_back();
		arr[1].pop_back();
		arr[2].pop_back();
		arr[3].pop_back();
	}
};