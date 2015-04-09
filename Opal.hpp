#pragma once

#include "Star_Constants.hpp"
#include "IFileHandle.hpp"
#include <vector>
#include <string>
#include "Graphable.hpp"
#include <iostream>
#include <sstream>
#include <exception>

class Star;

class OpalException : public exception{
public:

	OpalException(const char* str) :exception(str) {

	}

};

class Opal : public Graphable{
public:
	vector<vector<double>> table;
	Star* star;

	Opal(Star* star, string filename, bool use_opal);

	double get_kappa(double temperature, double density);

	void pushValues();

	inline virtual void popValues() {
		arr[0].pop_back();
		arr[1].pop_back();
	}
};