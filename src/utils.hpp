#pragma once


#define LOGERRORAT \
    "[ERROR] : " << __FILE__ << ":"  <<__LINE__ << " : "


#include <tuple>

#include "nanoflann/nanoflann.hpp"
#include "nanoflann/util.h"

#include "Color.hpp"    
#include "Vector2.hpp"
#include "ScatterPoint.hpp"
#include "DataSource.hpp"
#include "RandNumGenerator.hpp"
#include "reader.hpp"


namespace utils
{


int getInt(const char* argv);
DataSource<float> loadRawData(const char* volumePath);
void generateImage(std::string &filename, Image* image, int width, int height);
std::tuple<int, int> getBoundaryOfDataVavlue(float*& data, int dataCount);
std::tuple<Vector2<int>, Vector2<float>> convertImgIdxToVoxelIdx(int x, int y, const Vector2<float>& voxelImgRatio);
unsigned int calculateIdx(int x, int y, int z, int lateral);
bool isInsideSphere(const Point& center, float radius, const Point& p);
void generateSphereTestData(unsigned int lateral, float*& data);
const Point convertIdx1DTo3D(int idx, const int* dimension);
int convertIdx3DTo1D(const Point& p, const int* dimension);


template<typename DataType>
Point2<DataType> sampleDataRandom(const DataType* data, unsigned int dataCount, float ratio, ScatterPoint<DataType>*& sampledData)
{
    // std::cout << "[sampleDataRandom]\n";
    int numSample = static_cast<int>(dataCount*ratio);

    std::vector<int> sampleIndices;
    sampleIndices.reserve(dataCount);

    for(int i = 0; i < dataCount; i++)
    {
        sampleIndices.emplace_back(i);
    }

    std::shuffle(sampleIndices.begin(), sampleIndices.end(), std::default_random_engine(0));

    // std::cout << "[scatterData] (Before)\n";
    ScatterPoint<DataType>* scatterData = new ScatterPoint<DataType>[numSample];
    for(int i = 0; i < numSample; i++)
    {
        scatterData[i] = { sampleIndices[i], static_cast<DataType>(data[sampleIndices[i]]) };
    }
    sampledData = scatterData;
    // std::cout << "[scatterData] (After)\n";

    DataType maxValue = 0;
    DataType minValue = 0;
    DataType tmp = 0;

    for(int i = 0; i < numSample; i++)
    {
        tmp = scatterData[i].funcValue;
        if(maxValue < tmp)
        {
            maxValue = tmp;
        }
        if(minValue > tmp)
        {
            minValue = tmp;
        }
    }

    return { minValue, maxValue };
}


template<typename DataType>
void generatePointCloud(PointCloud<float>*& point, const DataSource<ScatterPoint<DataType>>& data)
{
    // std::cout << "[generatePointCloud]\n";
	point->pts.resize(data.count);

	for(unsigned int i = 0; i < data.count; i++)
	{
        Vector3 idx = convertIdx1DTo3D(data.data[i].index, data.dimension);

        point->pts[i].x = idx.x();
        point->pts[i].y = idx.y();
        point->pts[i].z = idx.z();
	}
    std::cout << std::endl;
}


template<typename Tree, typename DataType>
Tree* buildKDTree(const DataSource<ScatterPoint<DataType>>& data)
{
    PointCloud<float>* cloud = new PointCloud<float>;
    generatePointCloud(cloud, data);

    // construct a kd-tree index:
    // std::cout << "[kd-tree]\n";
    Tree* index = new Tree(3, *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index->buildIndex();

    return index;
}


}