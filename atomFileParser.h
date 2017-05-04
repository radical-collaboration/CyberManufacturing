#ifndef ATOMFILEPARSER_H
#define ATOMFILEPARSER_H


#include <vector>
#include <array>
#include <map>
//#include <pair>
#include <string>

typedef struct
{
    std::array <double, 3> velocity;
    std::vector <unsigned int> c_ccVec;
    double f_fpacc;
} collisionData;

//for collision files
typedef std::map<unsigned int, std::vector<collisionData>> mapCollisionData;
//for impact files
//first value for impact with wall
//second value for impact with impeller
typedef std::pair<unsigned int, unsigned int> pairImpactData;

mapCollisionData collisionFileParser (std::string filePath, std::string collisionFileName, double& time);

pairImpactData impactFileParser (std::string filePath, std::string impactFileName);
#endif // ATOMFILEDATA_H
