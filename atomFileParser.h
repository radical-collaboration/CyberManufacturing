#ifndef ATOMFILEPARSER_H
#define ATOMFILEPARSER_H


#include <vector>
#include <array>
#include <map>
#include <string>

typedef struct
{
    std::array <double, 3> velocity;
    //double c_cc;
    std::vector <unsigned int> c_ccVec;
    double f_fpacc;
} particleData;

typedef std::map<int, std::vector<particleData>> mapParticleIdData;

double getTimeValueFromAtomFileName (std::string fileName);

mapParticleIdData atomFileParser (std::string fileName);

#endif // ATOMFILEDATA_H
