#ifndef __PRECISION_H_RP__
#define __PRECISION_H_RP__

 #include <mkl_types.h>
 #include <mkl.h>
#include <mkl_vsl.h>
#include <float.h>
#include <math.h>
#ifndef SINGLE_PRECISION_RP__
#define vNumberRngUniform(...) vdRngUniform(__VA_ARGS__)
#define vNumberLn(...) vdLn(__VA_ARGS__)
#define cblas_numscal(...) cblas_dscal(__VA_ARGS__)
typedef double numeric_t_rp;
typedef double rate_t_rp;
typedef double time_t_rp;
typedef double effect_t_rp;
#define MEPS_RP	DBL_EPSILON
#define powfunc(...) pow(__VA_ARGS__)
#define vNumberRngUniform(...) vdRngUniform(__VA_ARGS__)
#define vNumberRngGaussian(...) vdRngGaussian(__VA_ARGS__)
#define vNumberRngExponential(...) vdRngExponential(__VA_ARGS__)
#define vNumberRngLognormal(...) vdRngLognormal(__VA_ARGS__)
#define vNumberRngGamma(...) vdRngGamma(__VA_ARGS__)
#define vNumberRngBeta(...) vdRngBeta(__VA_ARGS__)

#else
#define vNumberRngUniform(...) vsRngUniform(__VA_ARGS__)
#define vNumberLn(...) vsLn(__VA_ARGS__)
#define cblas_numscal(...) cblas_sscal(__VA_ARGS__)
typedef float rate_t_rp;
typedef float numeric_t_rp;
typedef float time_t_rp;
typedef float effect_t_rp;
#define powfunc(...) powf(__VA_ARGS__)
#define MEPS_RP	FLT_EPSILON
#define vNumberRngUniform(...) vsRngUniform(__VA_ARGS__)
#define vNumberRngGaussian(...) vsRngGaussian(__VA_ARGS__)
#define vNumberRngExponential(...) vsRngExponential(__VA_ARGS__)
#define vNumberRngLognormal(...) vsRngLognormal(__VA_ARGS__)
#define vNumberRngGamma(...) vsRngGamma(__VA_ARGS__)
#define vNumberRngBeta(...) vsRngBeta(__VA_ARGS__)
#endif

#endif