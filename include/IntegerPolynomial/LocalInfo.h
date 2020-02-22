#ifndef _LOCALINFO_H_
#define _LOCALINFO_H_


NTL_START_IMPL

typedef struct LocalInfo {
   long n;
   long NumPrimes;
   long NumFactors;
   vec_long p;
   vec_vec_long pattern;
   ZZ PossibleDegrees;
   PrimeSeq s;
} LocalInfoT;


NTL_END_IMPL

#endif
