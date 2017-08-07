#include <iostream>
#include <cmath>
#include <float.h>
#include <sstream>
#include <iomanip>
#include <dirent.h>
#include <cstring>
#include <omp.h>
#include "utility.h"

using namespace std;

void fun()
{
    cout << "utility::fun is called" << endl;
}

vectorOfInt2D getVectorOfInt2D(int n, int m, int val)
{
    vectorOfInt2D vectorInt2D(n, vector<int>(m, val));
    return vectorInt2D;
}

vectorOfInt3D getvectorOfInt3D(int n, int m, int p, int val)
{
    vectorOfInt3D vectorInt3D(n, vector<vector<int>>(m, vector<int>(p, val)));
    return vectorInt3D;
}

vectorOfInt4D getVectorOfInt4D(int n, int m, int p, int q, int val)
{
    vectorOfInt4D vectorInt4D(n, vector<vector<vector<int>>>(m, vector<vector<int>>(p, vector<int>(q, val))));
    return vectorInt4D;
}

vectorOfDouble2D getVectorOfDouble2D(int n, int m, double val)
{
    vectorOfDouble2D vectorDouble2D(n, vector<double>(m, val));
    return vectorDouble2D;
}

vectorOfDouble3D getVectorOfDouble3D(int n, int m, int p, double val)
{
    vectorOfDouble3D vectorDouble3D(n, vector<vector<double>>(m, vector<double>(p, val)));
    return vectorDouble3D;
}

vectorOfDouble4D getVectorOfDouble4D(int n, int m, int p, int q, double val)
{
    vectorOfDouble4D vectorDouble4D(n, vector<vector<vector<double>>>(m, vector<vector<double>>(p, vector<double>(q, val))));
    return vectorDouble4D;
}

double getMinimumOf3DVector(vectorOfDouble3D vector3D)
{
    double minValue = DBL_MAX;
    for (size_t d1 = 0; d1 < vector3D.size(); d1++)
        for (size_t d2 = 0; d2 < vector3D[d1].size(); d2++)
            for (size_t d3 = 0; d3 < vector3D[d1][d2].size(); d3++)
            {
                minValue = min(minValue, vector3D[d1][d2][d3]);
            }
    return minValue;
}

double getMinimumOf2DVector(vectorOfDouble2D vector2D)
{
    double minValue = DBL_MAX;
    for (size_t d1 = 0; d1 < vector2D.size(); d1++)
        for (size_t d2 = 0; d2 < vector2D[d1].size(); d2++)
        {
            minValue = min(minValue, vector2D[d1][d2]);
        }
    return minValue;
}

double getMaximumOfVector(std::vector<double> vec)
{
    double maxValue = -DBL_MAX;
    for (auto v : vec)
        maxValue = max(maxValue, v);

    return maxValue;
}

// double getMaximumOf2DVector(vectorOfDouble2D vector2D)
// {
//     double maxValue = -DBL_MAX;
//     for (size_t d1 = 0; d1 < vector2D.size(); d1++)
//         for (size_t d2 = 0; d2 < vector2D[d1].size(); d2++)
//         {
//             maxValue = max(maxValue, vector2D[d1][d2]);
//         }
//     return maxValue;
// }
double getMaximumOf2DArray(arrayOfDouble2D array2D)
{
    double maxValue = -DBL_MAX;
    for (size_t d1 = 0; d1 < array2D.size(); d1++)
        for (size_t d2 = 0; d2 < array2D[d1].size(); d2++)
        {
            maxValue = max(maxValue, array2D[d1][d2]);
        }
    return maxValue;
}

int getCountOfNegativeIn3DVector(vectorOfDouble3D vector3D)
{
    int count = 0;
    for (size_t d1 = 0; d1 < vector3D.size(); d1++)
        for (size_t d2 = 0; d2 < vector3D[d1].size(); d2++)
            for (size_t d3 = 0; d3 < vector3D[d1][d2].size(); d3++)
            {
                if (vector3D[d1][d2][d3] < 0.0)
                    count++;
            }
    return count;
}

vector<double> linearize3DVector(vectorOfDouble3D vector3D)
{
    // vector<double> data;
    size_t dim1 = vector3D.size();
    size_t dim2 = vector3D[0].size();
    size_t dim3 = vector3D[0][0].size();
    vector<double> data(dim1 * dim2 * dim3, 0.0);
    for (size_t d1 = 0; d1 < dim1; d1++)
        for (size_t d2 = 0; d2 < dim2; d2++)
            for (size_t d3 = 0; d3 < dim3; d3++)
                data[d1 * dim2 * dim3 + d2 * dim3 + d3] = vector3D[d1][d2][d3];

    return data;
}

string moreSigs(double d, int prec)
{
    std::stringstream ss;
    ss << d;
    ss.str("");
    ss << std::setprecision(prec) << std::fixed << d;

    std::string str;
    ss >> str;
    std::string::size_type s;
    for (s = str.length() - 1; s > 0; --s)
    {
        if (str[s] == '0')
            str.erase(s, 1);
        else
            break;
    }
    if (str[s] == '.')
        str.erase(s, 1);
    return str;
}

vector<string> listFiles(string path, string ext)
{
    bool gIgnoreHidden = true;
    string dotExt = "." + ext;
    vector<string> filelist;
    DIR *dirFile = opendir(path.c_str());
    if (dirFile)
    {
        struct dirent *hFile;
        errno = 0;
        while ((hFile = readdir(dirFile)) != NULL)
        {
            if (!strcmp(hFile->d_name, "."))
                continue;
            if (!strcmp(hFile->d_name, ".."))
                continue;

            // in linux hidden files all start with '.'
            if (gIgnoreHidden && (hFile->d_name[0] == '.'))
                continue;

            // dirFile.name is the name of the file. Do whatever string comparison
            // you want here. Something like:
            if (strstr(hFile->d_name, dotExt.c_str()))
                filelist.push_back(hFile->d_name);
            //cout<< "found an " << ext << " file " << hFile->d_name << endl;
        }
        closedir(dirFile);
    }
    return filelist;
}

arrayOfDouble2D vector2Array2D(vectorOfDouble2D vector2D)
{
   arrayOfDouble2D array2D;
   for (int i = 0; i < NUMBEROFFIRSTSOLIDBINS; i++)
   {
      vector<double> p = vector2D[i];
      std::copy (p.data(), p.data() + NUMBEROFSECONDSOLIDBINS, (array2D[i]).begin());
   }
   return array2D;
}

arrayOfInt2D vector2Array2D(vectorOfInt2D vector2D)
{
   
   arrayOfInt2D array2D;
   int i = 0;
   #pragma omp parallel default(shared) private(i)
   {
       #pragma omp for 
       for (i = 0; i < NUMBEROFFIRSTSOLIDBINS; i++)
       {
          vector<int> p = vector2D[i];
          std::copy (p.data(), p.data() + NUMBEROFSECONDSOLIDBINS, (array2D[i]).begin());
       }
   }
   return array2D;
}

vectorOfDouble2D array2Vector2D(arrayOfDouble2D array2D)
{
   vectorOfDouble2D vector2D = getVectorOfDouble2D(NUMBEROFFIRSTSOLIDBINS, NUMBEROFSECONDSOLIDBINS);
   int i = 0;
   #pragma omp parallel default(shared) private(i)
   {
      #pragma omp for 
      for (i = 0; i < NUMBEROFFIRSTSOLIDBINS; i++)
      std::copy ((array2D[i]).data(), (array2D[i]).data() + NUMBEROFSECONDSOLIDBINS, (vector2D[i]).begin());
    }
   return vector2D;

}

//arrayOfDouble3D vector2Array3D(vectorOfDouble3D vector3D);
//arrayOfDouble4D vector2Array4D(vectorOfDouble4D vector4D);
