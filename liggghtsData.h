#ifndef LIGGGHTSDATA_H
#define LIGGGHTSDATA_H

#include <map>
#include <vector>
#include "atomFileParser.h"
#include "utility.h"

typedef std::map <int, std::vector<particleData>> mapParticleIdData;
class liggghtsData
{
    static bool instanceFlag;
    static liggghtsData *lData;// = nullptr;
    liggghtsData();
    //void readDumpAtomFiles();

    std::map <double, mapParticleIdData> mapParticleIdDataOverTime;

public:
    //liggghtsData() = delete;
    static liggghtsData* getInstance();

    void readDumpAtomFiles();

    mapParticleIdData getMapParticledata (double time);
    arrayOfDouble2D getFinalNumberOfCollisions();
    ~liggghtsData();
};

#endif // LIGGGHTSDATA_H
