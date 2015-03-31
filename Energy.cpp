#include "Energy.hpp"
#include "Star.hpp"
#include <iostream>


long double Energy::get() {
	return _e_pp + _e_CNO + _e_3alpha;
}

 void Energy::pushValues(){
	 arr[0].push_back(get());
	 arr[1].push_back( _e_pp );
	 arr[2].push_back( _e_CNO );
	 arr[3].push_back( _e_3alpha );
	 
 }