#ifndef COMPARTMENT_H
#define COMPARTMENT_H

#include <vector>
#include "utility.h"
#include "parameters.h"

typedef struct
{
    vectorOfDouble2D fAll;
    vectorOfDouble2D fLiquid;
    vectorOfDouble2D fGas;
    double liquidAdditionRate;

    std::vector<double> vs;
    std::vector<double> vss;

    vectorOfDouble2D sMeshXY;
    vectorOfDouble2D ssMeshXY;

    vectorOfInt2D sAggregationCheck;
    vectorOfInt2D ssAggregationCheck;

    vectorOfInt2D sInd;
    vectorOfInt2D ssInd;

    vectorOfInt2D sIndB;
    vectorOfInt2D ssIndB;

    vectorOfDouble2D sLow;
    vectorOfDouble2D sHigh;

    vectorOfDouble2D ssLow;
    vectorOfDouble2D ssHigh;

    vectorOfInt2D sCheckB;
    vectorOfInt2D ssCheckB;

    vectorOfDouble2D diameter;

} CompartmentIn;

typedef struct
{
    vectorOfDouble2D dfAlldt;
    vectorOfDouble2D dfLiquiddt;
    vectorOfDouble2D dfGasdt;
    vectorOfDouble2D liquidBins;
    vectorOfDouble2D gasBins;
    vectorOfDouble2D internalVolumeBins;
    vectorOfDouble2D externalVolumeBins;
    vectorOfDouble4D aggregationKernel;
    vectorOfDouble4D breakageKernel;
} CompartmentOut;

typedef struct
{
    std::vector<double> DEMDiameter;
    vectorOfDouble2D DEMCollisionData;
    std::vector<double> DEMImpactData;
} CompartmentDEMIn;

typedef struct
{
    vectorOfDouble2D fAllPreviousCompartment;
    vectorOfDouble2D flPreviousCompartment;
    vectorOfDouble2D fgPreviousCompartment;
    vectorOfDouble2D fAllComingIn;
    vectorOfDouble2D fgComingIn;
} PreviousCompartmentIn;

CompartmentOut performCompartmentCalculations(PreviousCompartmentIn prevCompIn, CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double time, double timeStep);

#endif // COMPARTMENT_H
