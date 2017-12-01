#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <limits>

#define DOUBLEPRECISION std::numeric_limits<double>::epsilon()
#define FLOATPRECISION std::numeric_limits<float>::epsilon()
#define EPSILON 0.0 //DOUBLEPRECISION

//TUNING PARAMETERS
#define AGGREGATIONKERNELCONSTANT 1.0e-9//4.375e-7//3.75488281e-9//1.0//9.375e-4//1.0e-6 //1.0e-3 //240.0
//#define AGGREGATIONKERNELCONSTANTALPHA 0.9473
//#define AGGREGATIONKERNELCONSTANTBETA 0.5
//#define AGGREGATIONKERNELCONSTANTGAMMA 1.0
//#define AGGREGATIONKERNELCONSTANTDELTA 1.0
#define BREAKAGEKERNELCONSTANT 0.0 //5.0
#define BREAKAGEPROBABILITY 0.0
#define CONSOLIDATIONCONSTANT 0.0 //1.0e-5
//#define RPMCONSOLIDATIONCONSTANT 0.0
#define INITIALPOROSITY 0.5
#define MINIMUMPOROSITY 0.2
#define GRANULESATURATIONFACTOR 0.0

//MATERIAL & PROCESS PARAMETERS
#define GRANULATORLENGTH 0.38 //meter
#define NUMBEROFCOMPARTMENTS 4
#define DISTANCEBETWEENCOMPARTMENTS (GRANULATORLENGTH / NUMBEROFCOMPARTMENTS) // m
#define PARTICLERESIDENCETIME 49.0//490.0//20.0                                            //11.07 // seconds
#define PARTICLEAVERAGEVELOCITY (GRANULATORLENGTH / PARTICLERESIDENCETIME)

#define IMPELLERDIAMETER 0.114 // meter
#define IMPELLERSPEED 2000.0//1000.0   //2000.0 // RPM

#define PREMIXINGTIME 45.0//45.0                             //0.0 // seconds
#define LIQUIDADDITIONTIME 45.0//105.0//75.0//45.0
#define POSTMIXINGTIME 0.0//25.0                        // seconds
#define FINALTIME (PREMIXINGTIME + LIQUIDADDITIONTIME + POSTMIXINGTIME) // seconds
//#define TIMESTEP 1e-1 // seconds

#define THROUGHPUT 15.0    // kg / hr
#define SOLIDDENSITY 476.0 // kg / m^3
#define LIQUIDTOSOLIDRATIO 0.35
#define LIQUIDDENSITY 1000.0                                                              // kg / m^3
#define LIQUIDADDITIONRATE (((LIQUIDTOSOLIDRATIO * THROUGHPUT) / LIQUIDDENSITY) / 3600.0) // m3 / s
#define SHEARRATE (0.10472 * IMPELLERSPEED * (IMPELLERDIAMETER / 2))

// DEFINE BINS
#define NUMBEROFFIRSTSOLIDBINS 16 //Number of first solid bins
#define SCOEF 5.0e-16             // m^3
#define SBASE 3.0                 // Volume increment between two bins

#define NUMBEROFSECONDSOLIDBINS 16 //Number of second solid bins
#define SSCOEF 5.0e-16             // m^3
#define SSBASE 3.0                 // Volume increment between two bins

//DEM data
#define COLLISIONEFFICIENCYCONSTANT 0.01
#define TIMESTEPDEM 5.0e-7
#define NUMBEROFDEMBINS 16

#define DEMAGGREGATIONKERNELVALUE 0.0
#define DEMAGGREGATIONKERNELCONST 1.0

#define DEMBREAKAGEKERNELVALUE 0.0
#define DEMBREAKAGEKERNELCONST 1.0

#endif //PARAMETERS_H
