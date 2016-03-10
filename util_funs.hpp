/*
 * util_funs.hpp
 *
 *  Created on: Mar 5, 2016
 *      Author: yiminliu
 */

#ifndef UTIL_FUNS_HPP_
#define UTIL_FUNS_HPP_

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include "opca_bimodal.h"
#include <exception>
#include <stdexcept>
#include <boost/math/special_functions/erf.hpp>

using namespace::boost::math;

void TranformUniform2Normal(int l , vector<double> &m){
	for(int i = 0; i < l ; i++){
		m[i] = 1.4142135623 * erf_inv( 2*m[i] - 1 );
	}
}

void GeneratePerm(int Nc, double *x, double *perm){
	for(int i=0; i < Nc ; i++){
		perm[i] = exp(x[i]);
	}
}

void SaveData(string fn, int n, double* data){
	//cout << "Saving data to " << fn << endl;
	ofstream out;
	out.open(fn.c_str());

	if (out.is_open()){
		for (int i = 0; i < n; i++){
			out << data[i] << endl;
		}
		out.close();
	}else{
		throw runtime_error("Can not write file while saving data.");
	}
}

#endif /* UTIL_FUNS_HPP_ */
