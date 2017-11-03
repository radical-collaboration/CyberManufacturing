#include <cmath>

#include "kernel.h"
#include "utility.h"
#include "liggghtsData.h"
using namespace std;

#define DUMP(varName) dumpData(varName, #varName)
#define DUMP2D(varName) dump2DData(varName, #varName)
#define DUMP3D(varName) dump3DData(varName, #varName)

#define DUMPCSV(varName) dumpCSV(varName, #varName)
#define DUMP2DCSV(varName) dump2DCSV(varName, #varName)
#define DUMP4DCSV(varName) dump4DCSV(varName, #varName)


arrayOfDouble4D DEMDependentAggregationKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, arrayOfDouble2D externalLiquidContent, double timeStep)
{
    arrayOfDouble4D aggregationKernel = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Initialization
    arrayOfDouble4D collisionFrequency = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
    arrayOfDouble4D collisionEfficiency = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //Collision Frequency (from 2D Number of Collisions)
   
    arrayOfDouble2D fAll = compartmentIn.fAll;

    arrayOfDouble2D DEMCollisionData = compartmentDEMIn.DEMCollisionData;

    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    collisionFrequency[s1][ss1][s2][ss2] = (DEMCollisionData[s2][ss2] * timeStep) / TIMESTEPDEM;

    // Constant Collision Efficiency
    // FROM Barrasso, Tamrakar, Ramachandran. Procedia Engineering 102 (2015) 1295�1304. (p. 1298)
    //    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
    //        for (int ss1  = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
    //            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
    //                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
    //                    collisionEfficiency[s1][ss1][s2][ss2] = COLLISIONEFFICIENCYCONSTANT;
    // Calculating critical velocity below which aggregation will occur

    liggghtsData* lData = liggghtsData::getInstance();
    vector<double> diaCol = lData->getDEMParticleDiameters();

    if((diaCol).size() == 0)
    {
        cout << "Velocity data is missing in LIGGGHTS collision file" << endl;
        return aggregationKernel;
    }

    vector<double> velocity = lData->getFinalDEMCollisionVelocity();
    if ((velocity).size() == 0)
    {
        cout << "Velocity data is missing in LIGGGHTS collision file" << endl;
        return aggregationKernel;
    }
    double inverseDiameterSum = 0.0;
    double inverseMassSum = 0.0;
    int sized = diaCol.size();

    for (int i = 0; i < sized; i++)
    {
    	inverseDiameterSum += (1/diaCol[i]);
    	inverseMassSum += (1 / (4 / 3) * M_PI * pow((diaCol[i] / 2),3) * SOLIDDENSITY);
    }

    double sum = 0.0;
    int size1 = velocity.size();
    double harmonic_diameter = sized / inverseDiameterSum;
    double harmonic_mass = sized / inverseMassSum;
    double Ucritical = (1 + (1/COEFFICIENTOFRESTITUTION)) * log((LIQUIDTHICKNESS / SURFACEASPERITIES)) * (3 * M_PI * pow(harmonic_diameter, 2) * BINDERVISCOSITY) / (8 * harmonic_mass);
    for (int i = 0; i < size1; i++)
    	sum += velocity[i];

    double averageVelocity = sum / NUMBEROFDEMBINS;
    double stdDevVelocity = 0.0;
    double varianceVelocity = 0.0;
    
    for (int i = 0; i < size1; ++i)
    {
        varianceVelocity += pow((velocity[i] - averageVelocity), 2) / 10;
    }
    
    stdDevVelocity = sqrt(varianceVelocity);
    //double intVelocity = 0.0;
    vector<double> probablityOfVelocity(size1, 0.0);
    for (int i = 0; i < size1; i++)
    {
     	probablityOfVelocity[i] = (1 / (velocity[i] * sqrt(2 * M_PI) * stdDevVelocity)) * exp(-((log(velocity[i]) - averageVelocity) / (2 * pow(varianceVelocity, 2))));
     	// cout << "Probability at " << velocity[i] << "is " << probablityOfVelocity[i] << endl;
    }

    /*for (int i = 0; i < size1 - 1; ++i)
    {
    	if( (velocity[i]) < Ucritical)
     		intVelocity += ((probablityOfVelocity[i + 1] + probablityOfVelocity[i]) / 2) * (velocity[i + 1] - velocity[i]);
    }*/

    // FROM Sen, Barrasso, Singh, Ramachandran. Processes 2014, 2, 89-111. (p. 96)
    // double collisionEfficiencyConstant=0.01;
    double criticalExternalLiquid = 0.2;
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                {
                    bool flag1 = (fAll[s1][ss1] >= 0.0) && (fAll[s2][ss2] >= 0.0);
                    bool flag2 = (externalLiquidContent[s1][ss1] >= criticalExternalLiquid) && (externalLiquidContent[s2][ss2] >= criticalExternalLiquid);
                    bool flag3 = (velocity[ss2] < Ucritical);
                    if (flag1 && flag2 && flag3)
                        //collisionEfficiency[s1][ss1][s2][ss2] = COLLISIONEFFICIENCYCONSTANT;
                        collisionEfficiency[s1][ss1][s2][ss2] = probablityOfVelocity[ss2];
                }
    
    //Aggregation Kernel Calculation
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                    aggregationKernel[s1][ss1][s2][ss2] = AGGREGATIONKERNELCONSTANT * collisionFrequency[s1][ss1][s2][ss2] * collisionEfficiency[s1][ss1][s2][ss2];

    return aggregationKernel;
}

arrayOfDouble4D DEMDependentBreakageKernel(CompartmentIn compartmentIn, CompartmentDEMIn compartmentDEMIn, double timeStep)
{
    arrayOfDouble4D breakageKernel = getArrayOfDouble4D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS, NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    //INITIALIZATION
    //arrayOfDouble2D impactFrequency = getArrayOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);

    arrayOfDouble2D fAll = compartmentIn.fAll;
    vector<double> numberOfImpacts = compartmentDEMIn.DEMImpactData;
    vector<double> impactFrequency = compartmentDEMIn.DEMImpactData;   
    // double Ubreak = 0.0;
    //Impact Frequency (from 1D Number of Impacts)
    for (int s = 0; s < NUMBEROFFIRSTSOLIDBINS; s++)
        for (int ss = 0; ss < NUMBEROFSECONDSOLIDBINS; ss++)
            for (int i = 0; i < NUMBEROFDEMBINS; i++)
            {
                if (fAll[s][ss] > 0.0)
                    //impactFrequency[s][ss] = (numberOfImpacts[i] * timeStep) / TIMESTEPDEM;
                    impactFrequency[i] = (numberOfImpacts[i] * timeStep) / TIMESTEPDEM;
            }

    // Critical velocity for breakage
    //cout << "CRITICALSTOKESDEFNUMBER = " << (2 * CRITICALSTOKESDEFNUMBER / SOLIDDENSITY) << endl;
    //cout << "last fraction = " << (BINDERVISCOSITY / compartmentIn.diameter[0][0]) << endl;  
    //cout << "ratio of porosity = " << (pow((1 - INITIALPOROSITY),2) / pow(INITIALPOROSITY,2)) << endl;
    double Ubreak = (2 * CRITICALSTOKESDEFNUMBER / SOLIDDENSITY) * (9 / 8.0) * (pow((1 - INITIALPOROSITY),2) / pow(INITIALPOROSITY,2)) * (9 / 16.0) * (BINDERVISCOSITY / compartmentIn.diameter[0][0]); 
    //cout << "Ubreak = " << Ubreak << endl;
    liggghtsData* lData = liggghtsData::getInstance();
    vector<double> velocity = lData->getFinalDEMImpactVelocity();
    if ((velocity).size() == 0)
    {
        cout << "Velocity data is missing in LIGGGHTS impact file" << endl;
        return breakageKernel;
    }
    
    double sum = 0.0;
    int size1 = velocity.size();

    for (int i = 0; i < size1; i++)
    	sum += velocity[i];

    double averageVelocity = sum / NUMBEROFDEMBINS;
    double stdDevVelocity = 0.0;
    double varianceVelocity = 0.0;
    for (int i = 0; i < size1; ++i)
    {
        varianceVelocity += pow((velocity[i] - averageVelocity), 2) / 10;
    }
    
    stdDevVelocity = sqrt(varianceVelocity);
    //double intVelocity = 0.0;
    // cout << "Std Dev. of Velocity = " << stdDevVelocity << endl;

    vector<double> probablityOfVelocity(size1, 0.0);
    for (int i = 0; i < size1; i++)
    {
     	probablityOfVelocity[i] = (1 / (velocity[i] * sqrt(2 * M_PI) * stdDevVelocity)) * exp(-((log(velocity[i]) - averageVelocity) / (2 * pow(varianceVelocity, 2))));
     	// cout << "Probability at " << velocity[i] << "is " << probablityOfVelocity[i] << endl;
    }

    /*for (int i = 0; i < size1 - 1; ++i)
    {
    	if( (velocity[i] / 100.0) > Ubreak)
     		intVelocity += ((probablityOfVelocity[i + 1] + probablityOfVelocity[i]) / 2) * (velocity[i + 1] - velocity[i]);
    
    }*/

    //DUMP(impactFrequency);
    
    //Breakage Kernel Calculation
    for (int s1 = 0; s1 < NUMBEROFFIRSTSOLIDBINS; s1++)
        for (int ss1 = 0; ss1 < NUMBEROFSECONDSOLIDBINS; ss1++)
            for (int s2 = 0; s2 < NUMBEROFFIRSTSOLIDBINS; s2++)
                for (int ss2 = 0; ss2 < NUMBEROFSECONDSOLIDBINS; ss2++)
                	if(velocity[ss2] > Ubreak)
                		breakageKernel[s1][ss1][s2][ss2] = impactFrequency[ss1] * probablityOfVelocity[ss2] * BREAKAGEKERNELCONSTANT;
                    //breakageKernel[s1][ss1][s2][ss2] = impactFrequency[ss1] * BREAKAGEPROBABILITY * BREAKAGEKERNELCONSTANT;
                    //breakageKernel[s1][ss1][s2][ss2] = impactFrequency[s1][ss1] * BREAKAGEPROBABILITY;

    return breakageKernel;
    
}
