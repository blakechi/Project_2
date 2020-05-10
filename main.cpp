#include <iostream>
#include <string>

#include "nanoflann/nanoflann.hpp"
#include "nanoflann/util.h"
#include "Eigen/Dense"

#include "utils.hpp"
#include "Interpolation.hpp"
#include "Timer.hpp"
#include "ScatterPoint.hpp"
#include "DataSource.hpp"


typedef nanoflann::KDTreeSingleIndexAdaptor
        <
            nanoflann::L2_Simple_Adaptor<float, PointCloud<float>>,
            PointCloud<float>,
            3 /* dim */
		> KDTree;

typedef float DataType;


int main(int argc, char* argv[])
{

    std::string outputName = "image.ppm";
    std::string dataPath = argv[1];
    DataSource data = utils::loadRawData(dataPath.c_str());

    float sampleRate = static_cast<float>(utils::getInt(argv[2]))/100;
    ScatterPoint<DataType>* sampleData;
    utils::sampleDataRandom<DataType>(data.data, data.count, sampleRate, sampleData);
    DataSource<ScatterPoint<DataType>> scatterData;
    scatterData.data = sampleData;
    scatterData.dimension = data.dimension;
    scatterData.count = static_cast<int>(data.count*sampleRate);


    KDTree* kdTree = utils::buildKDTree<KDTree>(scatterData);


    Vector2<int> imageDimension = { 
        utils::getInt(argv[3]),
        utils::getInt(argv[4])
    };
    Vector2<float> voxelImgRatio = { 
        static_cast<float>(data.dimension[0])/imageDimension.x,
        static_cast<float>(data.dimension[1])/imageDimension.y
    };
    Image* image = new Color[imageDimension.x*imageDimension.y];


    std::cout << "[Main] Start to generate image" << std::endl;
    {
        Timer timer;

        for(int y = imageDimension.y - 1; y >= 0; y--)
        {
            for(int x = 0; x < imageDimension.x; x++)
            {
                
            }
        }
    }

    
    utils::generateImage(outputName, image, imageDimension.x, imageDimension.y);


    delete[] image;
    delete[] scatterData.data;
    delete[] scatterData.dimension;


    DataSource<ScatterPoint<DataType>> data_;
    int count = 27 - 1;
    ScatterPoint<DataType>* d = new ScatterPoint<DataType>[count];

    float tmp;
    for(int j = 0; j < count; j++)
    {
        if(j == 13)
        {
            continue;
        }
        else
        {
            tmp = static_cast<DataType>(j);
            d[j] = { j, tmp };
        }
    }
    data_.data = d;
    data_.count = count;
    data_.dimension = new int[3];
    data_.dimension[0] = 3;
    data_.dimension[1] = 3;
    data_.dimension[2] = 3;

    KDTree* kdTree_ = utils::buildKDTree<KDTree, DataType>(data_);
    Interpolation::localShepard2<KDTree>(kdTree_, data_, Vector3(1, 1, 1), 6);

    for(int j = 0; j < count; j++)
    {
        std::cout << data_.data[j].funcValue << " ";
    }
    std::cout << std::endl;

    delete[] data_.data;
    delete[] data_.dimension;
}