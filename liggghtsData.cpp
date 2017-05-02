#include <iostream>

#include "liggghtsData.h"
#include "parameters.h"


#define ATOMFILEPATH "./sampledumpfiles/"
using namespace std;

bool liggghtsData::instanceFlag = false;
liggghtsData* liggghtsData::lData = nullptr;

liggghtsData::liggghtsData()
{

}

liggghtsData* liggghtsData::getInstance()
{
    if(!instanceFlag)
    {
        //readDumpAtomFiles();
        lData = new liggghtsData();
        instanceFlag = true;
        return lData;
    }
    else
        return lData;
}

void liggghtsData::readDumpAtomFiles()
{
    if (!mapParticleIdDataOverTime.empty())
        return;

    vector <string> fileList = listFiles(ATOMFILEPATH, "atom");

    for(auto fileName : fileList)
    {
        //cout << fileName << endl;
        double time = getTimeValueFromAtomFileName(fileName);
        if (time == 0.0)
            continue;
        mapParticleIdData mapData = atomFileParser(ATOMFILEPATH + fileName);
        if (mapData.empty())
        {
            cout << fileName << " file is invalid" << endl;
            continue;
        }
        pair<double, mapParticleIdData> mapEntry(time, mapData);
        mapParticleIdDataOverTime.insert(mapEntry);
    }

//    cout << endl;
//    for(auto it = mapParticleIdDataOverTime.begin(); it != mapParticleIdDataOverTime.end(); it++)
//    {
//        //auto it = mapParticleIdDataOverTime.end();
//        //it--;
//        cout << it->first << endl;
//        for(auto it2 = (it->second).begin(); it2 != (it->second).end(); it2++)
//        {
//            cout << "particle type = " << it2->first << endl;
//            vector<particleData> pDataVec = it2->second;
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

}

mapParticleIdData liggghtsData::getMapParticledata (double time)
{
    mapParticleIdData mapPData;

    if (mapParticleIdDataOverTime.empty())
        return mapPData;

    auto it = mapParticleIdDataOverTime.find(time);
    if (it != mapParticleIdDataOverTime.end())
        mapPData = it->second;

    return mapPData;
}

arrayOfDouble2D liggghtsData::getFinalNumberOfCollisions()
{
    arrayOfDouble2D nCollisions;

    if(!instanceFlag)
        return nCollisions;

    if (mapParticleIdDataOverTime.empty())
        return nCollisions;

    auto mapIt = mapParticleIdDataOverTime.end();

    mapParticleIdData mapPData = getMapParticledata((--mapIt)->first);

    if (mapPData.empty())
        return nCollisions;

    nCollisions = getArrayOfDouble2D(NUMBEROFDEMBINS, NUMBEROFDEMBINS);

    for(auto itMapPData = mapPData.begin(); itMapPData != mapPData.end(); itMapPData++)
    {
        int row = itMapPData->first;

        vector<particleData> vecPData = itMapPData->second;

        for (auto pData : vecPData)
        {
            vector <unsigned int> c_ccVec = pData.c_ccVec;
            int c_ccCount = 0;
            for(auto c_cc : c_ccVec)
                nCollisions[row-1][c_ccCount++] += c_cc;
        }
    }
    return nCollisions;
}

