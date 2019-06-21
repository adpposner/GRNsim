//parameters_p.h - Private structures for model parameters 
// Only for network generation and for default values when desired
/*
    
    Copyright (C) 2018, Russell Posner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
#ifndef PARAMETERS_P_H_RP__
#define PARAMETERS_P_H_RP__
#include "globals.h"




struct ConnectionParametersUnit {
	int deg_low;
	int deg_high;
	effect_t_rp effect_prob;
	enum ConnDegreeTp byInOrOut;
	species_t connTp;
};

typedef struct connectionParameters {
	struct ConnectionParametersUnit *miRConns;
	struct ConnectionParametersUnit *tfCodingConns;
	struct ConnectionParametersUnit *tfNoncodingConns;
} connectionParameters;

 struct MessRates {
	rate_t_rp production;
	rate_t_rp decay;
} ;

 struct TFRates {
	rate_t_rp production;
	rate_t_rp decay;
	rate_t_rp affinity;
	rate_t_rp frequency;
	effect_t_rp effectMin;
	effect_t_rp effectMax;
} ;

 struct MicroRates {
	rate_t_rp production;
	rate_t_rp decay;
	rate_t_rp affinity;
	rate_t_rp frequency;
	effect_t_rp effectStrength;
} ;


struct defaultRates {
	struct MessRates mess;
	struct TFRates tf;
	struct MicroRates micro;
	rate_t_rp defaultRate;
};

struct initialQuantities {
	ulong_type coding;
	ulong_type noncoding;
	ulong_type mess;
	ulong_type micro;
	ulong_type TF;
};

struct ParameterStruct {
	connectionParameters connParms;
	struct defaultRates rates;
	struct initialQuantities initialQtys;
} system_parameters;


#endif
