#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "parameters.h"
#include "parameterData.h"
#define FILEPATH "./"
using namespace std;

bool parameterData::instanceFlag = false;
parameterData *parameterData::pData = nullptr;

parameterData::parameterData()
{
    cout << "LIGGGHTS data pointer created" << endl;
}

parameterData::~parameterData()
{
    cout << "LIGGGHTS data pointer destroyed" << endl;
    delete pData;
}

parameterData *parameterData::getInstance()
{
    if (!instanceFlag)
    {
        pData = new parameterData();
        instanceFlag = true;
        return pData;
    }
    else
        return pData;
}

void parameterData::readPBMInputFile (/*parameter to specify iteration specific file */)
{
    //For now, initialising values from parameters.h
    //Later, when format of input file is finalized, following values will be initialized by reading PBM input files

    aggKernelConst = AGGREGATIONKERNELCONSTANT; //Aggregation Kernel Constant
    brkKernelConst = BREAKAGEKERNELCONSTANT; //Breakage Kernel Constant
    consConst = CONSOLIDATIONCONSTANT; //Consolidation Constant
    initPorosity = INITIALPOROSITY; // Initial Porosity
    minPorosity = MINIMUMPOROSITY; //Minimum Porosity
    granSatFactor = GRANULESATURATIONFACTOR; //Granule Stauration Factor

    //Probablistic kernel constants
    critStDefNum = CRITICALSTOKESDEFNUMBER; //Critical Stokes Deformation number value
    bindVisc = BINDERVISCOSITY; //Binder viscosity
    coefOfRest = COEFFICIENTOFRESTITUTION; //Coefficient of restitution
    liqThick = LIQUIDTHICKNESS; //Liquid thickness
    surfAsp = SURFACEASPERITIES; //Surface asperities

    //Material & process parameters
    granulatorLength = GRANULATORLENGTH; //unit: meter
    nCompartments = NUMBEROFCOMPARTMENTS; //Number of compartmets
    partticleResTime = PARTICLERESIDENCETIME; //Particle Residence Time (unit: seconds)
    impDiameter = IMPELLERDIAMETER; //Diammeter of Impeller (unit: meter)
    impSpeed = IMPELLERSPEED; //Speed of impeller (unit: rpm)

    //Time Parameters (unit: seconds)
    premixTime = PREMIXINGTIME; //Pre mixing time
    liqAddTime = LIQUIDADDITIONTIME; //Liquid addition time
    postMixTime = POSTMIXINGTIME; //Post mixing time 

    //Other Process parameters
    throughput = THROUGHPUT; //unit: Kg/hr
    solDensity = SOLIDDENSITY; //Solid density (unit: Kg/m^3)
    liqSolidRatio = LIQUIDTOSOLIDRATIO; // Liquid to Solid ratio
    liqDensity = LIQUIDDENSITY; //Liquid density

    //Bins parameter
    //First particle
    nFirstSolidBins = NUMBEROFFIRSTSOLIDBINS;
    fsVolCoeff = SCOEF; //Volume of first solid bin
    fsVolBase = SBASE;  //Base value to increase volume between two bins
    //Second particle
    nSecondSolidBins = NUMBEROFSECONDSOLIDBINS;
    ssVolCoeff = SSCOEF; //Volume of second solid bin
    ssVolBase = SSBASE;  //Base value to increase volume between two bins

    //DEM data
    colEffConst = COLLISIONEFFICIENCYCONSTANT; // Collision efficiency constant
    demTimeStep = TIMESTEPDEM;
    nDEMBins = NUMBEROFDEMBINS;
}