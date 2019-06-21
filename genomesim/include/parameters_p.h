#ifndef PARAMETERS_P_H_RP__
#define PARAMETERS_P_H_RP__
#include "globals.h"
#include "netgen/networkdistributions.h"



struct ConnectionParametersUnit {
    CompoundDiscreteDistribution *connectionDistribution;
	enum ConnDegreeTp byInOrOut;
	species_t connTp;
};



typedef struct connectionParameters {
	struct ConnectionParametersUnit *miRConns;
	struct ConnectionParametersUnit *tfCodingConns;
	struct ConnectionParametersUnit *tfNoncodingConns;
} connectionParameters;



typedef struct ConstOrCtsDistribution{
	int isConst;
	union{
		rate_t_rp rate;
		CompoundContinuousDistribution* distn;
	};
} ConstOrCtsDistribution;



 struct MessRates {
	ConstOrCtsDistribution *production;
	ConstOrCtsDistribution *decay;
} ;

 struct TFRates {
	ConstOrCtsDistribution *production;
	ConstOrCtsDistribution *decay;
	ConstOrCtsDistribution *affinity;
	ConstOrCtsDistribution *frequency;
	ConstOrCtsDistribution *effectSize;
} ;

 struct MicroRates {
	ConstOrCtsDistribution *production;
	ConstOrCtsDistribution *decay;
	ConstOrCtsDistribution *affinity;
	ConstOrCtsDistribution *frequency;
	ConstOrCtsDistribution *prodEffectStrength;
    ConstOrCtsDistribution *decayEffectStrength;
} ;


struct defaultRates {
	struct MessRates mess;
	struct TFRates tf;
	struct MicroRates micro;
	ConstOrCtsDistribution * defaultRate;
};



struct initialQuantities {
	qty_type coding;
	qty_type noncoding;
	qty_type mess;
	qty_type micro;
	qty_type TF;
};

struct ParameterStruct {
	connectionParameters connParms;
	struct defaultRates rates;
	struct initialQuantities initialQtys;
};



#endif
