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


std::tuple<int, int> getBoundaryOfDataVavlue(int*& data, int dataCount)
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


DataSource<int> loadRawData(const char* volumePath)
{
    std::cout << "[Load] Reading data...\n";

    int   data_type, data_count;
    int*  data_dim = new int[3];
    void* data_ptr_void = NULL;
    ReadVolume(volumePath, data_type, data_count, data_dim[0], data_dim[1], data_dim[2], data_ptr_void);
    char* data_ptr_char = (char*)data_ptr_void;

    int* data_ptr = new int[data_count];
    for(int i = 0; i < data_count; i++)
    {
        data_ptr[i] = static_cast<int>(data_ptr_char[i]) + 128;
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


void generateSphereTestData(unsigned int lateral, int*& data)
{
    Point center(32, 32, 32);
    float radius = 20.0f;

    int* data_ = new int[lateral*lateral*lateral];

    for(int k  = 0; k < lateral; k++)
    {
        for(int j  = 0; j < lateral; j++)
        {
            for(int i  = 0; i < lateral; i++)
            {
                unsigned int idx = calculateIdx(i, j, k, lateral);

                if(isInsideSphere(center, radius, {i, j, k})) 
                { 
                    data_[idx] = 127; 
                }
                else
                {
                    data_[idx] = -128; 
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


void sampleDataRandom(const int* data, unsigned int dataCount, float ratio, ScatterPoint<int>*& sampledData)
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
    ScatterPoint<int>* scatterData = new ScatterPoint<int>[numSample];
    for(int i = 0; i < numSample; i++)
    {
        scatterData[i] = { sampleIndices[i], data[sampleIndices[i]] };
    }
    sampledData = scatterData;
    // std::cout << "[scatterData] (After)\n";
}

const Vector3 convertIdx1DTo3D(int idx, const int* dimension)
{
    // std::cout << "[convertIdx1DTo3D]\n";
    return Vector3(
        (idx % (dimension[0]*dimension[1])) % dimension[0],
        (idx % (dimension[0]*dimension[1]))/dimension[0],
        idx/(dimension[0]*dimension[1])
    );
}

void generatePointCloud(PointCloud<float> &point, const DataSource<ScatterPoint<int>>& data)
{
    // std::cout << "[generatePointCloud]\n";
	point.pts.resize(data.count);

	for(unsigned int i = 0; i < data.count; i++)
	{
        Vector3 idx = convertIdx1DTo3D(data.data[i].index, data.dimension);

        point.pts[i].x = idx.x();
        point.pts[i].y = idx.y();
        point.pts[i].z = idx.z();
	}
    std::cout << std::endl;
}


void interpolateData(const DataSource<ScatterPoint<int>>& data, const Point& target, const int numNeighbors)
{
    // std::cout << "[interpolateData]\n";
    PointCloud<float> cloud;
    generatePointCloud(cloud, data);

    // construct a kd-tree index:
    // std::cout << "[kd-tree]\n";
	typedef nanoflann::KDTreeSingleIndexAdaptor<
		nanoflann::L2_Simple_Adaptor<float, PointCloud<float>>,
		PointCloud<float>,
		3 /* dim */
		> kd_tree;

    kd_tree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
	index.buildIndex();


    float queryPt[3] = { 
        target.x(), 
        target.y(), 
        target.z() 
    };

    size_t numResults = numNeighbors;
    std::vector<size_t> nearestNeighbors(numResults);
    std::vector<float> neighborsDistance(numResults);

    // std::cout << "[knnSearch]\n";
    numResults = index.knnSearch(&queryPt[0], numResults, &nearestNeighbors[0], &neighborsDistance[0]);


    nearestNeighbors.resize(numResults);
    neighborsDistance.resize(numResults);
    
    // TODO: interpolation
    for(size_t neighbor : nearestNeighbors)
    {
        std::cout << "(" << neighbor << ", ";
        std::cout << data.data[neighbor].index << ", ";
        std::cout << data.data[neighbor].funcValue << ") ";
    }
    std::cout << std::endl;

    for(size_t d : neighborsDistance)
    {
        std::cout << d << " ";
    }
    std::cout << std::endl;
}


}