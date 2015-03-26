#pragma once

#include <algorithm>
#include "Graphable.hpp"

class Star;

class Energy : public Graphable{
public:
	Star* star;

	long double get();

	Energy(Star* star) : star(star), Graphable(1) {}

	virtual void pushValues() {
		arr[0].push_back(get());
	}
};