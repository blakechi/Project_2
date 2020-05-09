#include <iostream>
#include <string>

#include "utils.hpp"
#include "Interpolation.hpp"
#include "Timer.hpp"


int main(int argc, char* argv[])
{
    std::string outputName = "image.ppm";
    std::string dataPath = argv[1];
    utils::DataSource data = utils::loadRawData(dataPath.c_str());

    float sampleRate = static_cast<float>(utils::getInt(argv[2]))/100;
    utils::ScatterPoint<int>* sampleData;
    utils::sampleDataRandom(data.data, data.count, sampleRate, sampleData);
    utils::DataSource<utils::ScatterPoint<int>> scatterData;
    scatterData.data = sampleData;
    scatterData.dimension = data.dimension;
    scatterData.count = static_cast<int>(data.count*sampleRate);

    
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
                int gray = 0;
                for(int z = 0; z < data.dimension[2]; z++)
                {
                    gray += data.data[x + y*data.dimension[0] + z*data.dimension[0]*data.dimension[1]];
                }
                gray /= (data.dimension[2]);
                image[x + y*imageDimension.x] = Color(gray, gray, gray);
            }
        }
    }

    
    utils::generateImage(outputName, image, imageDimension.x, imageDimension.y);


    delete[] image;
    delete[] data.data;
    delete[] data.dimension;


    utils::DataSource<utils::ScatterPoint<int>> data_;
    int count = 27 - 1;
    utils::ScatterPoint<int>* d = new utils::ScatterPoint<int>[count];
    for(int j = 0; j < count; j++)
    {
        if(j == 13)
        {
            continue;
        }
        else
        {
            d[j] = { j, j };
        }
    }
    data_.data = d;
    data_.count = count;
    data_.dimension = new int[3];
    data_.dimension[0] = 3;
    data_.dimension[1] = 3;
    data_.dimension[2] = 3;

    utils::interpolateData(data_, Vector3(1, 1, 1), 6);

    for(int j = 0; j < count; j++)
    {
        std::cout << data_.data[j].funcValue << " ";
    }
    std::cout << std::endl;

    delete[] data_.data;
    delete[] data_.dimension;
}