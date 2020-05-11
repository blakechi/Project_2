#include <fstream>
#include <cmath>
#include <random>
#include <algorithm>
#include <string>
#include <stdexcept>


#include "utils.hpp"


#define BLANK -999


namespace utils
{


int getInt(const char* argv)
{
    std::string arg = argv;
    try 
    {
        std::size_t pos;
        int x = std::stoi(arg, &pos);
        if (pos < arg.size()) 
        {
            std::cerr << LOGERRORAT << "Trailing characters after number: " << arg << '\n';
        }

        return x;
    } 
    catch (std::invalid_argument const &ex) 
    {
        std::cerr << LOGERRORAT << "Invalid number: " << arg << '\n';
    } 
    catch (std::out_of_range const &ex) 
    {
        std::cerr << LOGERRORAT << "Number out of range: " << arg << '\n';
    }

    return -1;
}


void generateImage(std::string &fileName, Image* image, int width, int height)
{
    std::ofstream output(fileName, std::ios::binary);

    output << "P3\n" << width << " " << height << "\n255\n";

    for(int g = height - 1; g >= 0; g--)
    {
        for(int h = 0; h < width; h++)
        {
            output << image[g*width + h] << "\n";
        }
    }

    output.clear();
    output.close();
}


std::tuple<int, int> getBoundaryOfDataVavlue(float*& data, int dataCount)
{
    int maxValue = 0;
    int minValue = 0;
    int tmp = 0;

    for(int i = 0; i < dataCount; i++)
    {
        tmp = data[i];
        if(maxValue < tmp)
        {
            maxValue = tmp;
        }
        if(minValue > tmp)
        {
            minValue = tmp;
        }
    }

    return { maxValue, minValue };
}


DataSource<float> loadRawData(const char* volumePath)
{
    std::cout << "[Load] Reading data...\n";

    int   data_type, data_count;
    int*  data_dim = new int[3];
    void* data_ptr_void = NULL;
    ReadVolume(volumePath, data_type, data_count, data_dim[0], data_dim[1], data_dim[2], data_ptr_void);
    char* data_ptr_char = (char*)data_ptr_void;

    float* data_ptr = new float[data_count];
    for(int i = 0; i < data_count; i++)
    {
        data_ptr[i] = static_cast<float>(data_ptr_char[i]) + 128;
    }
    delete[] data_ptr_char;

    auto[data_max_value, data_min_value] = getBoundaryOfDataVavlue(data_ptr, data_count);
    std::cout << "[Load] Successfully read data\n";

    return { data_ptr, data_dim, data_count, data_max_value, data_min_value };
}


std::tuple<Vector2<int>, Vector2<float>> convertImgIdxToVoxelIdx(int x, int y, const Vector2<float>& voxelImgRatio)
{
    Vector2<int> index = { static_cast<int>(x*voxelImgRatio.x), static_cast<int>(y*voxelImgRatio.y) };
    Vector2<float> indexInVoxel = { (x*voxelImgRatio.x) - index.x, (y*voxelImgRatio.y) - index.y };
    return { index, indexInVoxel };
}


unsigned int calculateIdx(int x, int y, int z, int lateral)
{
    return x + y*lateral + z*lateral*lateral;
}


bool isInsideSphere(const Point& center, float radius, const Point& p)
{
    return ((center - p).magnitude() <= radius)? true : false;
}


void generateSphereTestData(unsigned int lateral, float*& data)
{
    Point center(lateral/2, lateral/2, lateral/2);
    float radius = 0.6*lateral;

    float* data_ = new float[lateral*lateral*lateral];

    for(int k  = 0; k < lateral; k++)
    {
        for(int j  = 0; j < lateral; j++)
        {
            for(int i  = 0; i < lateral; i++)
            {
                unsigned int idx = calculateIdx(i, j, k, lateral);

                if(isInsideSphere(center, radius, {i, j, k})) 
                { 
                    data_[idx] = 1 - ((center - Point(i, j, k)).magnitude()/radius); 
                }
                else
                {
                    data_[idx] = 0; 
                }
            }
        }
    }
    data = data_;

    // Example:
    // unsigned int lateral = 64;
    // char* data;
    // utils::generateSphereTestData(lateral, data);
    // int* dataDimension = new int[3];
    // dataDimension[0] = lateral;
    // dataDimension[1] = lateral;
    // dataDimension[2] = lateral;
}


const Point convertIdx1DTo3D(int idx, const int* dimension)
{
    // std::cout << "[convertIdx1DTo3D]\n";
    return Point(
        (idx % (dimension[0]*dimension[1])) % dimension[0],
        (idx % (dimension[0]*dimension[1]))/dimension[0],
        idx/(dimension[0]*dimension[1])
    );
}


int convertIdx3DTo1D(const Point& p, const int* dimension)
{
    return p.x() + p.y()*dimension[0] + p.z()*dimension[0]*dimension[1];
}


}