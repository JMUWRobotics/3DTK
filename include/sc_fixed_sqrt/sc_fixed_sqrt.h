#ifndef SC_FIXED_SQRT_H
#define SC_FIXED_SQRT_H

#define SC_INCLUDE_FX

#define WL 16
#define IWL 8

#include <systemc.h>


using namespace sc_dt;
using namespace std;

typedef sc_fixed<WL, IWL> fixed_t;

fixed_t sc_fixed_heron_sqrt(fixed_t s);


#endif // SC_FIXED_SQRT_H