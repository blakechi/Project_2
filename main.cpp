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
    DataSource<float> data = utils::loadRawData(dataPath.c_str());

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


    int lateral = 32;
    DataSource<float> data_;

    utils::generateSphereTestData(lateral, data_.data);
    data_.count = lateral*lateral*lateral;
    data_.dimension = new int[3];
    data_.dimension[0] = lateral;
    data_.dimension[1] = lateral;
    data_.dimension[2] = lateral;

    ScatterPoint<DataType>* sampleData_;
    utils::sampleDataRandom<DataType>(data_.data, data_.count, sampleRate, sampleData_);

    DataSource<ScatterPoint<DataType>> data__;
    data__.count = static_cast<int>(data_.count*sampleRate);
    data__.dimension = data_.dimension;
    data__.data = sampleData_;

    KDTree* kdTree_ = utils::buildKDTree<KDTree, DataType>(data__);


    Image* image_ = new Color[imageDimension.x*imageDimension.y];
    for(int i = 0; i < imageDimension.x*imageDimension.y; i++)
    {
        image_[i] = Color(0);
    }

    // random sampled
    for(int idx = 0; idx < data__.count; idx++)
    {
        Point p = utils::convertIdx1DTo3D(data__.data[idx].index, data__.dimension);
        if(p.z() == 16)
        {
            image_[static_cast<int>(p.x() + p.y()*data__.dimension[0])] = Color(data__.data[idx].funcValue).convert255();
        }
    }

    std::string imageName = "sampled.ppm";
    utils::generateImage(imageName, image_, imageDimension.x, imageDimension.y);

    // interpolated
    for(int y = imageDimension.y - 1; y >= 0; y--)
    {
        for(int x = 0; x < imageDimension.x; x++)
        {
            // ScatterPoint<DataType> result = Interpolation::globalShepard2<KDTree>(kdTree_, data__, Vector3(x, y, 32), 6);
            ScatterPoint<DataType> result = Interpolation::localHardy<KDTree>(kdTree_, data__, Vector3(x, y, 32), 6, false);
            image_[x + y*imageDimension.x] = Color(result.funcValue).clamp().gammaCorrection(4).convert255();
        }
    }

    imageName = "interpolated.ppm";
    utils::generateImage(imageName, image_, imageDimension.x, imageDimension.y);

    delete[] data__.data;
    delete[] data__.dimension;
}