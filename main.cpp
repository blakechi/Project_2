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
#include "ColorMap.hpp"


#define NUMNEIGHBORS 6


float R = 10.0f;


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
    Point2<DataType> minMaxValue = utils::sampleDataRandom<DataType>(data.data, data.count, sampleRate, sampleData);
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
    for(int i = 0; i < imageDimension.x*imageDimension.y; i++)
    {
        image[i] = Color(0);
    }

    Image* image_ = new Color[imageDimension.x*imageDimension.y];


    // TODO: extand to multiple dimension
    int sliceIdx = utils::getInt(argv[5]);

    Color colorDummy(0);
    // Eigen::VectorXd CHardy = Interpolation::precomputeGlobalC<DataType>(data__, R);


    // original
    std::string imageName = "original.ppm";
    for(int y = imageDimension.y - 1; y >= 0; y--)
    {
        for(int x = 0; x < imageDimension.x; x++)
        {
            colorDummy = ColorMap::mapColorFrom(data.data[x + y*data.dimension[0] + sliceIdx*data.dimension[0]*data.dimension[1]], minMaxValue);
            image_[x + y*imageDimension.x] = colorDummy.clamp().convert255();
        }
    }
    utils::generateImage(imageName, image_, imageDimension.x, imageDimension.y);


    // random sampled
    imageName = "sampled.ppm";
    for(int idx = 0; idx < scatterData.count; idx++)
    {
        Point p = utils::convertIdx1DTo3D(scatterData.data[idx].index, scatterData.dimension);
        if(p.z() == sliceIdx)
        {
            colorDummy = ColorMap::mapColorFrom(scatterData.data[idx].funcValue, minMaxValue);
            image[static_cast<int>(p.x() + p.y()*imageDimension.x)] = colorDummy.clamp().convert255();
        }
    }
    utils::generateImage(imageName, image, imageDimension.x, imageDimension.y);


    // interpolated
    imageName = "interpolated.ppm";
    Timer timer(imageName);
    {
        float tmp = 0.0f;
        for(int y = imageDimension.y - 1; y >= 0; y--)
        {
            for(int x = 0; x < imageDimension.x; x++)
            {
                // ScatterPoint<DataType> result = Interpolation::localShepard2<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx), NUMNEIGHBORS);
                ScatterPoint<DataType> result = Interpolation::globalShepard2<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx));
                // ScatterPoint<DataType> result = Interpolation::localHardy<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx), false, R, NUMNEIGHBORS);
                // ScatterPoint<DataType> result = Interpolation::localHardy<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx), true, R, NUMNEIGHBORS);
                // ScatterPoint<DataType> result = Interpolation::approximateGlobalHardy<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx), false, R);
                // ScatterPoint<DataType> result = Interpolation::approximateGlobalHardy<KDTree, DataType>(kdTree, scatterData, Vector3(x, y, sliceIdx), true, R);
                // ScatterPoint<DataType> result = Interpolation::globalHardy<KDTree>(kdTree, CHardy, scatterData, Vector3(x, y, sliceIdx), false, R);
                // ScatterPoint<DataType> result = Interpolation::globalHardy<KDTree>(kdTree, CHardy, scatterData, Vector3(x, y, sliceIdx), true, R);
                image[x + y*imageDimension.x] = Color(result.funcValue).clamp();

                tmp = result.funcValue;
                if(minMaxValue.max < tmp)
                {
                    minMaxValue.max = tmp;
                }
                if(minMaxValue.min > tmp)
                {
                    minMaxValue.min = tmp;
                }

                if(x == 0 && y % 16 == 0)
                {
                    std::cout << y << '\n';
                }
            }
        }
    }

    for(int y = imageDimension.y - 1; y >= 0; y--)
    {
        for(int x = 0; x < imageDimension.x; x++)
        {
            colorDummy = ColorMap::mapColorFrom(image[x + y*imageDimension.x].r(), minMaxValue);
            image[x + y*imageDimension.x] = colorDummy.clamp().convert255();        
        }
    }       
    utils::generateImage(imageName, image, imageDimension.x, imageDimension.y);


    std::cout << "Averge Pixel Error: " << utils::calculateAveragePixelErrorBetween(image, image_, imageDimension.x*imageDimension.y) << '\n';


    delete[] image;
    delete[] scatterData.data;
    delete[] scatterData.dimension;
}