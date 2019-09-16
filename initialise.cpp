#include "Prot.h"
double dt[Nrpath];
double dt_minus_one[Nrpath];
void initalise_dt()
{
//This should be called once at the start of the program to initialise the global dt and dt_minus_one array
    for (int i = 0; i < (int) (Nrpath/2); ++i)
    {
        dt[i] = deltaT;
        dt[i + (int) (Nrpath/2)] = -1.0*deltaT;
    }
    
    dt_minus_one[0] = -1.0*deltaT;
    for (int i = 1; i < (int) (Nrpath/2 + 1); ++i)
    {
      dt_minus_one[i] = deltaT;
    }
    for (int i = (int) (Nrpath/2 + 1); i < Nrpath; ++i)
    {
      dt_minus_one[i] = -1.0*deltaT;
    }
}
