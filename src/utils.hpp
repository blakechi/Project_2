#pragma once


#define LOGERRORAT \
    "[ERROR] : " << __FILE__ << ":"  <<__LINE__ << " : "


#include <tuple>

#include "nanoflann/nanoflann.hpp"
#include "nanoflann/util.h"

#include "Color.hpp"    
#include "Vector2.hpp"
#include "RandNumGenerator.hpp"
#include "reader.hpp"


namespace utils
{


template<typename T>
struct ScatterPoint
{
    int index;
    T funcValue;
};


template<typename T>
struct DataSource
{
    T* data;
    int* dimension; // [x, y, z]
    int count;
    int maxValue;
    int minValue;
};

int getInt(const char* argv);
DataSource<int> loadRawData(const char* volumePath);
void generateImage(std::string &filename, Image* image, int width, int height);
std::tuple<int, int> getBoundaryOfDataVavlue(int*& data, int dataCount);
std::tuple<Vector2<int>, Vector2<float>> convertImgIdxToVoxelIdx(int x, int y, const Vector2<float>& voxelImgRatio);
unsigned int calculateIdx(int x, int y, int z, int lateral);
bool isInsideSphere(const Point& center, float radius, const Point& p);
void generateSphereTestData(unsigned int lateral, int*& data);
void sampleDataRandom(const int* data, unsigned int dataCount, float ratio, ScatterPoint<int>*& sampledData);
const Vector3 convertIdx1DTo3D(int idx, const int* dimension);
void generatePointCloud(PointCloud<float> &point, const DataSource<ScatterPoint<int>>& data);
void interpolateData(const DataSource<ScatterPoint<int>>& data, const Point& target, const int numNeighbors);

}