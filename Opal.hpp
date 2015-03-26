#pragma once

#include "Star_Constants.hpp"
#include <IFileHandle.hpp>
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

	Opal(Star* star, string filename) : star(star) {
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

	double get_kappa(double temperature, double density);

	void iterate();
};