#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include "atomFileParser.h"
using namespace std;


double getTimeValueFromAtomFileName (string fileName)
{
    const char * str = fileName.c_str();
    char digits[] = "1234567890";
    size_t firstDigitPos = strcspn (str,digits);
    if (firstDigitPos == fileName.length())
    {
        cout << fileName << " file name doesn't contain any time value" << endl;
        return 0.0;
    }
    size_t dotPos = fileName.find(".");
    if(dotPos == static_cast<size_t>(string::npos))
    {
        cout << fileName << " doesn't have any '.' for file extension" << endl;
        return 0.0;
    }
    string timeStr = fileName.substr (firstDigitPos, dotPos - firstDigitPos);
    return stod(timeStr);
}

mapParticleIdData atomFileParser (string fileName)
{
    mapParticleIdData mapIdData;
    ifstream atomFile;
    atomFile.open (fileName.c_str(), ifstream::in);
    if(!atomFile.is_open())
    {
        std::cout << "Unable to open " << fileName << " file" << endl;
        return mapIdData;
    }
    string line;
    string tmpStr;
    stringstream lineData;
    while (tmpStr.compare("ATOMS") && !atomFile.eof())
    {
        getline(atomFile, line);
        lineData = move(stringstream(line));
        lineData >> tmpStr;
        lineData >> tmpStr;
    }

    if(tmpStr.compare("ATOMS") || atomFile.eof())
    {
        //cout << fileName << " doesn't contain require info" << endl;
        return mapIdData;
    }

    while(tmpStr.compare("fz"))
        lineData >> tmpStr;

    if(tmpStr.compare("fz"))// || atomFile.eof())
    {
        //cout << fileName << " doesn't contain require info" << endl;
        return mapIdData;
    }

    int c_ccCount = 0;
    while(tmpStr.compare("f_fppacc"))
    {
        lineData >> tmpStr;
        c_ccCount++;
    }
    c_ccCount--;

    if(tmpStr.compare("f_fppacc"))// || atomFile.eof())
    {
        //cout << fileName << " doesn't contain require info" << endl;
        return mapIdData;
    }

    while (getline(atomFile, line))
    {
        lineData = move(stringstream(line));
        lineData >> tmpStr; //read and ignore particle id
        lineData >> tmpStr; //read particle type
        int particleType = stoi(tmpStr);
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore x, y & z value;
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore ix, iy & iz value;
        particleData pData;
        lineData >> pData.velocity[0] >> pData.velocity[1] >> pData.velocity[2]; //read vx, vy & vz value;
        lineData >> tmpStr >> tmpStr >> tmpStr; //read and ignore fx, fy & fz value;
        //read collision data
        pData.c_ccVec.resize(c_ccCount);
        for(int i = 0; i < c_ccCount; i++)
            lineData >> pData.c_ccVec[i];

        lineData >> pData.f_fpacc;

        auto mapIt = mapIdData.find(particleType);
        if( mapIt == mapIdData.end())
        {
            vector<particleData> tmpVec;
            tmpVec.push_back(pData);
            pair<int, vector<particleData>> mapEntry(particleType, tmpVec);
            mapIdData.insert(mapEntry);
        }
        else
            (mapIt->second).push_back(pData);

    }

//    cout << fileName << endl;
//    cout << mapIdData.size() << endl;
//    for(auto it = mapIdData.begin(); it != mapIdData.end(); it++)
//        cout << (it->second).size() << endl;
//    cout << endl;

    atomFile.close();
    return mapIdData;
}
