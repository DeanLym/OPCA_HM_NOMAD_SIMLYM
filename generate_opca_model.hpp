/*
 * generate_opca_model.hpp
 *
 *  Created on: Mar 5, 2016
 *      Author: yiminliu
 */

#ifndef GENERATE_OPCA_MODEL_HPP_
#define GENERATE_OPCA_MODEL_HPP_


#include "opca_bimodal.h"
#include <vector>
#include <exception>
#include <stdexcept>
using namespace std;


OPCA_BIMODAL* GenerateOPCAModel(){
	int l = 1, Nc = 1, Nr = 1;
	double gamma=1.0, mu1=1.0, mu2=1.0, v1=1.0, v2=1.0, logk_min=0.0, logk_max=1.0;
	string usig_file, xm_file;
	string temp;
	ifstream in;
	in.open("OPCA_PARAM.DATA");
	if(in.is_open()){
		in >> temp >> l;
		in >> temp >> Nc;
		in >> temp >> Nr;
		in >> temp >> gamma;
		in >> temp >> mu1;
		in >> temp >> mu2;
		in >> temp >> v1;
		in >> temp >> v2;
		in >> temp >> logk_min;
		in >> temp >> logk_max;
		in >> temp >> usig_file;
		in >> temp >> xm_file;
	}else{
		throw runtime_error("Can not open OPCA parameter file...");
	}
	in.close();
	OPCA_BIMODAL* opca_bm = new OPCA_BIMODAL(l,gamma,Nc,Nr); // declare an opca_binary object
	opca_bm->InputXm(xm_file.c_str());
	opca_bm->InputUSig(usig_file.c_str());
	opca_bm->set_bounds( logk_min , logk_max );
	opca_bm->set_mu( mu1 , mu2 );
	opca_bm->set_var( v1 , v2 );

	return opca_bm;
}


#endif /* GENERATE_OPCA_MODEL_HPP_ */
