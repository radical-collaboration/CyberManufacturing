#include "kernel.h"
#include "utility.h"
#include <omp.h>
#include <iostream>
using namespace std;

#define DUMP(varName) dumpData(varName, #varName)
#define DUMP2D(varName) dump2DData(varName, #varName)
#define DUMP3D(varName) dump3DData(varName, #varName)

#define DUMPCSV(varName) dumpCSV(varName, #varName)
#define DUMP2DCSV(varName) dump2DCSV(varName, #varName)
#define DUMP4DCSV(varName) dump4DCSV(varName, #varName)

arrayOfDouble4D DEMDependentAggregationKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, arrayOfDouble2D externalLiquidContent, double timeStep)
{
    arrayOfDouble4D aggregationKernel{{0.0}};//} = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Initialization
    //arrayOfDouble4D collisionFrequency{{0.0}};// = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
    //arrayOfDouble4D collisionEfficiency{{0.0}};// = getVectorOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Collision Frequency (from 2D Number of Collisions)
   
    arrayOfDouble2D fAll = vector2Array2D(compartmentIn.fAll);
    double criticalExternalLiquid = 0.2;
    int s1, ss1, s2, ss2;

    vectorOfDouble2D DEMCollisionData = compartmentDEMIn.DEMCollisionData;
    //arrayOfDouble2D DEMCollisionData = vector2Array2D(compartmentDEMIn.DEMCollisionData);
    #pragma omp parallel for default(shared) private(s1, ss1, s2, ss2) 
    for (s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                {
                    //collisionFrequency[s1][ss1][s2][ss2] = (DEMCollisionData[s2][ss2] * timeStep) / TIMESTEPDEM;
                    double collisionFrequency = (DEMCollisionData[s2][ss2] * timeStep) / TIMESTEPDEM;
                    double collisionEfficiency = 0.0;
                    bool flag1 = (fAll[s1][ss1] >= 0.0) && (fAll[s2][ss2] >= 0.0);
                    bool flag2 = (externalLiquidContent[s1][ss1] >= criticalExternalLiquid) && (externalLiquidContent[s2][ss2] >= criticalExternalLiquid);
                    if (flag1 && flag2)                        
                        collisionEfficiency/*[s1][ss1][s2][ss2]*/ = COLLISIONEFFICIENCYCONSTANT;

                    #pragma omp atomic write
                    aggregationKernel[s1][ss1][s2][ss2] = AGGREGATIONKERNELCONSTANT * collisionFrequency/*[s1][ss1][s2][ss2]*/ * collisionEfficiency/*[s1][ss1][s2][ss2]*/;
                }
    return aggregationKernel;
}

arrayOfDouble4D DEMDependentBreakageKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double timeStep)
{
    arrayOfDouble4D breakageKernel{{0.0}};

    //INITIALIZATION
    arrayOfDouble2D impactFrequency{{0.0}};
    vectorOfDouble2D fAll = compartmentIn.fAll;
    vector<double> numberOfImpacts = compartmentDEMIn.DEMImpactData;
   
    //Impact Frequency (from 1D Number of Impacts)
    for (int s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
        for (int ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            for (int i = 0; i < NUMBEROFDEMBINS - 1; i++)
            {
                if (fAll[s][ss] > 0.0)
                    impactFrequency[s][ss] = (numberOfImpacts[i] * timeStep) / (fAll[s][ss] * TIMESTEPDEM);
            }

    //Breakage Kernel Calculation
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    breakageKernel[s1][ss1][s2][ss2] = impactFrequency[s1][ss1] * BREAKAGEPROBABILITY;

    return breakageKernel;
}
