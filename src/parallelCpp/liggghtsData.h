#ifndef LIGGGHTSDATA_H
#define LIGGGHTSDATA_H

#include <map>
#include <vector>
//#include <pair>
#include "atomFileParser.h"
#include "utility.h"

class liggghtsData
{
    static bool instanceFlag;
    static liggghtsData *lData; // = nullptr;
    liggghtsData();
    bool checkFileConsistency(std::string collisionFile, std::string impactFile);

    std::map<double, mapCollisionData> mapCollisionDataOverTime;
    std::map<double, pairImpactData> mapImpactDataOverTime;
    //std::map <double, mapParticleDiameter> mapParticleDiameterOverTime;

  public:
    //liggghtsData() = delete;
    static liggghtsData *getInstance();

    void readLiggghtsDataFiles();

    mapCollisionData getMapCollisionData(double time);
    pairImpactData getPairImpactData(double time);

    arrayOfDouble2D getFinalNumberOfCollisions();
    std::vector<double> getFinalNumberOfImpacts();
    std::vector<double> getParticleDiameters();
    ~liggghtsData();
};

#endif // LIGGGHTSDATA_H
