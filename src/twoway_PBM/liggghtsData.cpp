#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>

#include "liggghtsData.h"
#include "parameterData.h"

#define ATOMFILEPATH "./sampledumpfiles/"
using namespace std;

bool liggghtsData::instanceFlag = false;
liggghtsData *liggghtsData::lData = nullptr;

liggghtsData::liggghtsData()
{
    cout << "LIGGGHTS data pointer created" << endl;
}

liggghtsData::~liggghtsData()
{
    cout << "LIGGGHTS data pointer destroyed" << endl;
    delete lData;
}

liggghtsData *liggghtsData::getInstance()
{
    if (!instanceFlag)
    {
        //readDumpAtomFiles();
        lData = new liggghtsData();
        instanceFlag = true;
        return lData;
    }
    else
        return lData;
}

bool liggghtsData::checkFileConsistency(std::string collisionFile, std::string impactFile)
{
    const char *collisionFileStr = collisionFile.c_str();
    char digits[] = "1234567890";
    size_t firstDigitPos = strcspn(collisionFileStr, digits);
    if (firstDigitPos == collisionFile.length())
    {
        cout << collisionFile << " file name doesn't contain any time value" << endl;
        return false;
    }
    size_t dotPos = collisionFile.find(".");
    if (dotPos == static_cast<size_t>(string::npos))
    {
        cout << collisionFile << " doesn't have any '.' for file extension" << endl;
        return false;
    }
    string timeStr = collisionFile.substr(firstDigitPos, dotPos - firstDigitPos);
    double timeInCollisionFile = abs(stod(timeStr));

    const char *impactFileStr = impactFile.c_str();
    firstDigitPos = strcspn(impactFileStr, digits);
    if (firstDigitPos == impactFile.length())
    {
        cout << impactFile << " file name doesn't contain any time value" << endl;
        return false;
    }
    dotPos = impactFile.find(".");
    if (dotPos == static_cast<size_t>(string::npos))
    {
        cout << impactFile << " doesn't have any '.' for file extension" << endl;
        return false;
    }
    timeStr = impactFile.substr(firstDigitPos, dotPos - firstDigitPos);
    double timeInImpactFile = stod(timeStr);

    return (timeInCollisionFile == timeInImpactFile);
}

void liggghtsData::readLiggghtsDataFiles(string coreVal, string diaVal)
{
    if (!mapCollisionDataOverTime.empty() && !mapImpactDataOverTime.empty())
        return;

    // vector<string> fileList = listFiles(ATOMFILEPATH, "atom");
    string fileExt = coreVal + string("_") + diaVal;
    vector<string> fileList = listFiles(ATOMFILEPATH, fileExt);

    string subStrCollision = "collision";
    string subStrImpact = "impact";

    size_t lengthSubStrCollision = subStrCollision.size();
    size_t lengthSubStrImpact = subStrImpact.size();

    vector<string> collisionFileList;
    vector<string> impactFileList;

    for (auto fileName : fileList)
    {
        if (!(fileName.substr(0, lengthSubStrCollision)).compare(subStrCollision))
            collisionFileList.push_back(fileName);

        if (!(fileName.substr(0, lengthSubStrImpact)).compare(subStrImpact))
            impactFileList.push_back(fileName);
    }

    size_t nCountCollisionFiles = collisionFileList.size();
    size_t nCountImpactFiles = impactFileList.size();

    if (nCountCollisionFiles != nCountImpactFiles)
    {
        cout << "Collision & Impact files are not in sync" << endl;
        return;
    }

    sort(collisionFileList.begin(), collisionFileList.end());
    sort(impactFileList.begin(), impactFileList.end());

    fileList.clear();

    for (size_t c = 0; c < nCountCollisionFiles; c++)
    {
        auto collisionFile = collisionFileList[c];
        auto impactFile = impactFileList[c];

        bool check = checkFileConsistency(collisionFile, impactFile);
        if (!check)
            continue;

        double time = 0.0;
        mapParticleIdToType mapPartIdToType;
        mapCollisionData mapColData = collisionFileParser(ATOMFILEPATH, collisionFile, time, mapPartIdToType);
        if (mapColData.empty())
        {
            cout << collisionFile << " file is invalid" << endl;
            continue;
        }

        mapImpactData mapImpData = impactFileParser(ATOMFILEPATH, impactFile, mapPartIdToType);
        if (mapImpData.empty())
        {
            cout << impactFile << " file is invalid" << endl;
            continue;
        }

        pair<double, mapCollisionData> mapCollisionEntry(time, mapColData);
        mapCollisionDataOverTime.insert(mapCollisionEntry);

        pair<double, mapImpactData> mapImpactEntry(time, mapImpData);
        mapImpactDataOverTime.insert(mapImpactEntry);
    }

    //    cout << endl;
    //    for(auto it = mapCollisionDataOverTime.begin(); it != mapCollisionDataOverTime.end(); it++)
    //    {
    //        cout <<"timestep = " <<  it->first << endl;
    //        for(auto it2 = (it->second).begin(); it2 != (it->second).end(); it2++)
    //        {
    //            cout << "particle type = " << it2->first << endl;
    //            cout << "particle diameter = " << get<0>(it2->second) << endl;
    //            vector<collisionData> pDataVec = get<1>(it2->second);
    //            cout << "no. of rows = " << pDataVec.size() << endl;
    //            for (auto vec : pDataVec)
    //            {
    //                cout << "vx = " << vec.velocity[0] << endl;
    //                cout << "vy = " << vec.velocity[1] << endl;
    //                cout << "vz = " << vec.velocity[2] << endl;
    //                for (size_t i = 1; i <= (vec.c_ccVec).size(); i++)
    //                    cout << "c_cc_" << i << '\t';
    //                cout << endl;
    //                for(auto c_cc : vec.c_ccVec)
    //                    cout << c_cc << '\t';
    //                cout << endl;
    //                cout << "f_fpacc = " << vec.f_fpacc << endl << endl;
    //            }
    //        }
    //
    //    }
    //
    //    cout << endl;
    //    for(auto it = mapImpactDataOverTime.begin(); it != mapImpactDataOverTime.end(); it++)
    //    {
    //        cout <<"timestep = " <<  it->first << endl;
    //        cout << "Number of impact with wall = " << (it->second).first << endl;
    //        cout << "Number of impact with impeller = " << (it->second).second << endl << endl;
    //    }
}

mapCollisionData liggghtsData::getMapCollisionData(double time)
{
    mapCollisionData mapData;

    if (mapCollisionDataOverTime.empty())
        return mapData;

    auto it = mapCollisionDataOverTime.find(time);
    if (it != mapCollisionDataOverTime.end())
        mapData = it->second;

    return mapData;
}

mapImpactData liggghtsData::getMapImpactData(double time)
{
    mapImpactData mapData;
    if (mapImpactDataOverTime.empty())
       return mapData;

    auto it = mapImpactDataOverTime.find(time);
    if (it != mapImpactDataOverTime.end())
        mapData = it->second;

    return mapData;
}

vector<double> liggghtsData::getFinalDEMImpactData()
{
    vector<double> nImpacts;

    if (!instanceFlag)
        return nImpacts;

    if (mapImpactDataOverTime.empty())
        return nImpacts;

    auto mapIt = mapImpactDataOverTime.end();

    mapImpactData mapData = getMapImpactData((--mapIt)->first);

    if (mapData.empty())
    return nImpacts;
    
    parameterData *pData = parameterData::getInstance();
    unsigned int nDEMBins = pData->nDEMBins;
    
    nImpacts.resize(nDEMBins);

    for (auto itMapData = mapData.begin(); itMapData != mapData.end(); itMapData++)
        {
            int row = itMapData->first;
            nImpacts[row-1] = (itMapData->second).size();
        }

    return nImpacts;
}

arrayOfDouble2D liggghtsData::getFinalDEMCollisionData()
{
    arrayOfDouble2D nCollisions;

    if (!instanceFlag)
        return nCollisions;

    if (mapCollisionDataOverTime.empty())
        return nCollisions;

    auto mapIt = mapCollisionDataOverTime.end();

    mapCollisionData mapData = getMapCollisionData((--mapIt)->first);

    if (mapData.empty())
        return nCollisions;

    parameterData *pData = parameterData::getInstance();
    unsigned int nDEMBins = pData->nDEMBins;
    
    nCollisions = getArrayOfDouble2D(nDEMBins, nDEMBins);
    
    vector<size_t> particleTypeCount;
    for (auto itMapData = mapData.begin(); itMapData != mapData.end(); itMapData++)
    {
        int row = itMapData->first;

        vector<collisionData> vecCollisionData = get<1>(itMapData->second);

        particleTypeCount.push_back(vecCollisionData.size());
        for (auto data : vecCollisionData)
        {
            vector<int> c_ccVec = data.c_ccVec;
            int c_ccCount = 0;
            for (auto c_cc : c_ccVec)
                nCollisions[row - 1][c_ccCount++] += c_cc;
        }
    }

    for (size_t i = 0; i < nCollisions.size(); i++)
        for (size_t j = 0; j < nCollisions[i].size(); j++)
        {
            nCollisions[i][j] = nCollisions[i][j] / (particleTypeCount[i] * particleTypeCount[j]);
            nCollisions[i][j] = std::isnan(nCollisions[i][j]) ? 0.0 : nCollisions[i][j];
        }

    return nCollisions;
}

vector<double> liggghtsData::getDEMParticleDiameters()
{
    vector<double> particleDiameters;

    if (!instanceFlag)
        return particleDiameters;

    if (mapCollisionDataOverTime.empty())
        return particleDiameters;

    auto mapIt = mapCollisionDataOverTime.end();

    mapCollisionData mapData = getMapCollisionData((--mapIt)->first);

    if (mapData.empty())
        return particleDiameters;

    for (auto itMapData = mapData.begin(); itMapData != mapData.end(); itMapData++)
        particleDiameters.push_back(get<0>(itMapData->second));

    return particleDiameters;
}

vector<double> liggghtsData::getFinalDEMCollisionVelocity()
{
	vector<double> velocityCollision;
	vector<double> velocityIntCollision;

	if (!instanceFlag)
        return velocityCollision;

    if (mapCollisionDataOverTime.empty())
        return velocityCollision;

    auto mapIt = mapCollisionDataOverTime.end();

    for(int i = 0; i < 3; i++)
    {
	    mapCollisionData mapData = getMapCollisionData((--mapIt)->first);

	    if(mapData.empty())
	    	continue;
        
        parameterData *pData = parameterData::getInstance();
        unsigned int nDEMBins = pData->nDEMBins;
	    
        velocityCollision.resize(nDEMBins);
	    velocityIntCollision.resize(nDEMBins);

	    for(auto itMapData = mapData.begin(); itMapData != mapData.end(); itMapData++)
	    {
            array<double,3> aveColVel{{0.0}};
	    	vector<collisionData> vecColData = get<1>(itMapData->second);
            size_t nPartilcesOfEachType = vecColData.size();
	    	int row = itMapData->first;
	    	for(auto vecData : vecColData)
	    	{
	    		aveColVel[0] += fabs(vecData.velocity[0]);
	    		aveColVel[1] += fabs(vecData.velocity[1]);
	    		aveColVel[2] += fabs(vecData.velocity[2]);
	    	}
	    	velocityIntCollision[row - 1] = sqrt(pow(aveColVel[0],2) + pow(aveColVel[1],2) + pow(aveColVel[2],2)) / nPartilcesOfEachType;
	    	velocityCollision[row -1] += velocityIntCollision[row -1] / 3;
	    }
	}
    return velocityCollision;

}

vector<double> liggghtsData::getFinalDEMImpactVelocity()
{
    vector<double> velocity;
    vector<double> velocityInt;
    if (!instanceFlag)
        return velocity;

    if (mapImpactDataOverTime.empty())
        return velocity;

    auto mapIt = mapImpactDataOverTime.end();
    

    for (int i = 0; i < 3; ++i)
    {
	    mapImpactData mapData = getMapImpactData((--mapIt)->first);

	    if (mapData.empty())
	    //return velocity;
	    	continue;
        array<double, 3> aveVeloComp{{0.0}};

        parameterData *pData = parameterData::getInstance();
        unsigned int nDEMBins = pData->nDEMBins;
	    
        velocity.resize(nDEMBins);
	    velocityInt.resize(nDEMBins);

	    for (auto itMapData = mapData.begin(); itMapData != mapData.end(); itMapData++)
	    {
            vector<impactData> vecImpData = itMapData->second;
            size_t nPartilcesOfEachType = vecImpData.size();
            int row = itMapData->first;            
            for (auto impData : vecImpData)
            {
                aveVeloComp[0] += fabs(impData.velocity[3]);
                aveVeloComp[1] += fabs(impData.velocity[4]);
                aveVeloComp[2] += fabs(impData.velocity[5]);
            }
            velocityInt[row - 1] = sqrt(pow(aveVeloComp[0], 2) + pow(aveVeloComp[1], 2) + pow(aveVeloComp[2], 2)) / nPartilcesOfEachType / 1000;
            velocity[row -1] += velocityInt[row -1] / 3; 
        }
    }

    //velocity[row - 1] = sqrt(pow(aveVeloComp[0], 2) + pow(aveVeloComp[1], 2) + pow(aveVeloComp[2], 2));
    return velocity;
}