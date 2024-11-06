#ifndef _STG_INPUTS_H
#define _STG_INPUTS_H
#include<string>

struct  STGInputs{
   static double delta0;
   static double ub0;
   static double lx, ly, lz;
   static double xinflow, yinflow, zinflow;
   static bool lGenerateSEM;
   static double xlinf;
   static double qinfd, rhoinf;
   static void read_input();
};

#endif // _STG_INPUTS_H


