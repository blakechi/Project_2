#pragma once


#include "Eigen/Dense"

#include "utils.hpp"
#include "Vector3.hpp"
#include "DataSource.hpp"
#include "ScatterPoint.hpp"


namespace Interpolation
{


// Tri-linear
//
//
template<typename T>
float triLinear3D(const T* funcValue, const Point& p)
{
    float x = p.x();
    float y = p.y();
    float z = p.z();

    return (1 - x)*(1 - y)*(1 - z)*funcValue[0] // 000
        + (1 - x)*(1 - y)*z*funcValue[1] // 001
        + (1 - x)*y*(1 - z)*funcValue[2] // 010
        + (1 - x)*y*z*funcValue[3] // 011
        + x*(1 - y)*(1 - z)*funcValue[4] // 100
        + x*(1 - y)*z*funcValue[5] // 101
        + x*y*(1 - z)*funcValue[6] // 110
        + x*y*z*funcValue[7]; // 111
}

template<>
float triLinear3D<char>(const char* funcValue, const Point& p)
{
    float x = p.x();
    float y = p.y();
    float z = p.z();

    return (1 - x)*(1 - y)*(1 - z)*static_cast<int>(funcValue[0])
        + (1 - x)*(1 - y)*z*static_cast<int>(funcValue[1])
        + (1 - x)*y*(1 - z)*static_cast<int>(funcValue[2])
        + (1 - x)*y*z*static_cast<int>(funcValue[3])
        + x*(1 - y)*(1 - z)*static_cast<int>(funcValue[4])
        + x*(1 - y)*z*static_cast<int>(funcValue[5])
        + x*y*(1 - z)*static_cast<int>(funcValue[6])
        + x*y*z*static_cast<int>(funcValue[7]);
}

template<typename T>
Vector3 dTriLinear3D(const T* funcValue, const Point& p)
{
    float x = p.x();
    float y = p.y();
    float z = p.z();

    return Vector3(
        (-1)*(1 - y)*(1 - z)*funcValue[0]
        + (-1)*(1 - y)*z*funcValue[1]
        + (-1)*y*(1 - z)*funcValue[2]
        + (-1)*y*z*funcValue[3]
        + (1 - y)*(1 - z)*funcValue[4]
        + (1 - y)*z*funcValue[5]
        + y*(1 - z)*funcValue[6]
        + y*z*funcValue[7],
        
        (1 - x)*(-1)*(1 - z)*funcValue[0]
        + (1 - x)*(1 - z)*funcValue[1]
        + (1 - x)*(1 - z)*funcValue[2]
        + (1 - x)*z*funcValue[3]
        + x*(-1)*(1 - z)*funcValue[4]
        + x*(-1)*z*funcValue[5]
        + x*(-1)*z*funcValue[6]
        + x*z*funcValue[7],

        (1 - x)*(1 - y)*(-1)*funcValue[0]
        + (1 - x)*(1 - y)*funcValue[1]
        + (1 - x)*y*(-1)*funcValue[2]
        + (1 - x)*y*funcValue[3]
        + x*(1 - y)*(-1)*funcValue[4]
        + x*(1 - y)*funcValue[5]
        + x*y*(-1)*funcValue[6]
        + x*y*funcValue[7]
    );
}

template<>
Vector3 dTriLinear3D<char>(const char* funcValue, const Point& p)
{
    float x = p.x();
    float y = p.y();
    float z = p.z();

    return Vector3(
        (-1)*(1 - y)*(1 - z)*static_cast<int>(funcValue[0])
        + (-1)*(1 - y)*z*static_cast<int>(funcValue[1])
        + (-1)*y*(1 - z)*static_cast<int>(funcValue[2])
        + (-1)*y*z*static_cast<int>(funcValue[3])
        + (1 - y)*(1 - z)*static_cast<int>(funcValue[4])
        + (1 - y)*z*static_cast<int>(funcValue[5])
        + y*(1 - z)*static_cast<int>(funcValue[6])
        + y*z*static_cast<int>(funcValue[7]),
        
        (1 - x)*(-1)*(1 - z)*static_cast<int>(funcValue[0])
        + (1 - x)*(1 - z)*static_cast<int>(funcValue[1])
        + (1 - x)*(1 - z)*static_cast<int>(funcValue[2])
        + (1 - x)*z*static_cast<int>(funcValue[3])
        + x*(-1)*(1 - z)*static_cast<int>(funcValue[4])
        + x*(-1)*z*static_cast<int>(funcValue[5])
        + x*(-1)*z*static_cast<int>(funcValue[6])
        + x*z*static_cast<int>(funcValue[7]),

        (1 - x)*(1 - y)*(-1)*static_cast<int>(funcValue[0])
        + (1 - x)*(1 - y)*static_cast<int>(funcValue[1])
        + (1 - x)*y*(-1)*static_cast<int>(funcValue[2])
        + (1 - x)*y*static_cast<int>(funcValue[3])
        + x*(1 - y)*(-1)*static_cast<int>(funcValue[4])
        + x*(1 - y)*static_cast<int>(funcValue[5])
        + x*y*(-1)*static_cast<int>(funcValue[6])
        + x*y*static_cast<int>(funcValue[7])
    );
}


// Tri-cubic
//
//
template<typename T>
float triCubic3D(const T* funcValue, const Point& p)
{
    float x = p.x();
    float y = p.y();
    float z = p.z();

    // TODO: Control Points[64]

    return (1 - x)*(1 - y)*(1 - z)*funcValue[0]
        + (1 - x)*(1 - y)*z*funcValue[1]
        + (1 - x)*y*(1 - z)*funcValue[2]
        + (1 - x)*y*z*funcValue[3]
        + x*(1 - y)*(1 - z)*funcValue[4]
        + x*(1 - y)*z*funcValue[5]
        + x*y*(1 - z)*funcValue[6]
        + x*y*z*funcValue[7];
}


// Shepard
//
//
template<typename Tree>
ScatterPoint<float> localShepard2(const Tree* kDTree, const DataSource<ScatterPoint<float>>& sampleData, const Point& target, const int numNeighbors)
{
    // std::cout << "[interpolateData]\n";
    int originalIdx = utils::convertIdx3DTo1D(target, sampleData.dimension);
    float queryPt[3] = { 
        target.x(), 
        target.y(), 
        target.z() 
    };

    size_t numResults = 1;
    std::vector<size_t> kNearestNeighbors(numResults);
    std::vector<float> neighborsDistance(numResults);

    // std::cout << "[knnSearch]\n";
    numResults = kDTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);

    if(numResults == 1 && sampleData.data[kNearestNeighbors[0]].index == originalIdx)
    {
        return {originalIdx, sampleData.data[kNearestNeighbors[0]].funcValue};
    }
    else
    {
        numResults = numNeighbors;
        kNearestNeighbors.resize(numResults);
        neighborsDistance.resize(numResults);

        // std::cout << "[knnSearch]\n";
        numResults = kDTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);
        kNearestNeighbors.resize(numResults);
        neighborsDistance.resize(numResults);
        

        // TODO: interpolation
        for(size_t neighbor : kNearestNeighbors)
        {
            std::cout << neighbor << " ";
        }
        std::cout << std::endl;

        for(float distance : neighborsDistance)
        {
            std::cout << distance << " ";
        }
        std::cout << std::endl;

        float numerator = 0;
        float denominator = 0;
        for(int i = 0; i < numResults; i++)
        {
            float distanceInv = 1/neighborsDistance[i];
            numerator += (sampleData.data[kNearestNeighbors[i]].funcValue)*distanceInv;
            denominator += distanceInv;
        }

        return { 
            originalIdx,
            numerator/denominator
        };
    }
}


template<typename Tree>
ScatterPoint<float> globalShepard2(const Tree* kdTree, const DataSource<ScatterPoint<float>>& sampleData, const Point& target, const int numNeighbors)
{
    // std::cout << "[interpolateData]\n";
    int originalIdx = utils::convertIdx3DTo1D(target, sampleData.dimension);
    float queryPt[3] = { 
        target.x(), 
        target.y(), 
        target.z() 
    };

    size_t numResults = 1;
    std::vector<size_t> kNearestNeighbors(numResults);
    std::vector<float> neighborsDistance(numResults);

    // std::cout << "[knnSearch]\n";
    numResults = kdTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);

    if(numResults == 1 && sampleData.data[kNearestNeighbors[0]].index == originalIdx)
    {
        return {originalIdx, sampleData.data[kNearestNeighbors[0]].funcValue};
    }
    else
    {
        float numerator = 0;
        float denominator = 0;
        for(int i = 0; i < sampleData.count; i++)
        {
            Point dataPoint = utils::convertIdx1DTo3D(sampleData.data[i].index, sampleData.dimension);
            float distanceInv = 1/(target - dataPoint).magnitude();
            numerator += (sampleData.data[i].funcValue)*distanceInv;
            denominator += distanceInv;
        }

        return { 
            originalIdx,
            numerator/denominator
        };
    }
}


// Hardy
//
//
template<typename Tree>
ScatterPoint<float> localHardy(Tree* kdTree, const DataSource<ScatterPoint<float>>& sampleData, const Point& target, const int numNeighbors, const bool INVERSE)
{
    // std::cout << "[interpolateData]\n";
    int originalIdx = utils::convertIdx3DTo1D(target, sampleData.dimension);
    float queryPt[3] = { 
        target.x(), 
        target.y(), 
        target.z() 
    };

    size_t numResults = 1;
    std::vector<size_t> kNearestNeighbors(numResults);
    std::vector<float> neighborsDistance(numResults);

    // std::cout << "[knnSearch]\n";
    numResults = kDTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);

    if(numResults == 1 && sampleData.data[kNearestNeighbors[0]].index == originalIdx)
    {
        return {originalIdx, sampleData.data[kNearestNeighbors[0]].funcValue};
    }
    else
    {
        int R = 10;
        Eigen::MatrixXd M;
        Eigen::MatrixXd C;
        Eigen::MatrixXd F;

        M.resize(numNeighbors, numNeighbors);
        C.resize(numNeighbors);
        F.resize(numNeighbors);
        
        numResults = numNeighbors;
        kNearestNeighbors.resize(numResults);
        neighborsDistance.resize(numResults);

        // std::cout << "[knnSearch]\n";
        numResults = kDTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);
        kNearestNeighbors.resize(numResults);
        neighborsDistance.resize(numResults);
        
        float distance;
        for(int i = 0; i < numResults; i++)
        {
            F(i) = sampleData.data[i].funcValue;

            Point rowNeighbor = utils::convertIdx1DTo3D(
                sampleData.data[kNearestNeighbors(i)].index,
                sampleData.dimension
            );

            for(int j = 0; j < numResults; j++)
            {
                distance = (rowNeighbor - utils::convertIdx1DTo3D(
                    sampleData.data[kNearestNeighbors(j)].index,
                    sampleData.dimension
                )).magnitude();

                M(i, j) = std::sqrtf(R*R + distance);

                if(INVERSE)
                {
                    M(i, j) = 1/M(i, j);
                }
            }
        }


        C = M.householderQr().solve(F);


        float interpolatedValue = 0;
        for(int i = 0; i < numResults; i++)
        {
            distance = (target - utils::convertIdx1DTo3D(
                sampleData.data[kNearestNeighbors(i)].index,
                sampleData.dimension
            )).magnitude();

            if(INVERSE)
            {
                interpolatedValue += C(i)/std::sqrtf(R*R + distance);
            }
            else
            {
                interpolatedValue += C(i)*std::sqrtf(R*R + distance);
            }
        }

        return { 
            originalIdx,
            interpolatedValue
        };
    }
}


template<typename Tree>
ScatterPoint<float> globalHardy(Tree* kdTree, const Eigen::MatrixXd& C, const DataSource<ScatterPoint<float>>& sampleData, const Point& target, const bool INVERSE)
{
    // std::cout << "[interpolateData]\n";
    int originalIdx = utils::convertIdx3DTo1D(target, sampleData.dimension);
    float queryPt[3] = { 
        target.x(), 
        target.y(), 
        target.z() 
    };

    size_t numResults = 1;
    std::vector<size_t> kNearestNeighbors(numResults);
    std::vector<float> neighborsDistance(numResults);

    // std::cout << "[knnSearch]\n";
    numResults = kDTree->knnSearch(&queryPt[0], numResults, &kNearestNeighbors[0], &neighborsDistance[0]);

    if(numResults == 1 && sampleData.data[kNearestNeighbors[0]].index == originalIdx)
    {
        return {originalIdx, sampleData.data[kNearestNeighbors[0]].funcValue};
    }
    else
    {
        int R = 10;
        float distance;
        float interpolatedValue = 0;
        for(int i = 0; i < sampleData.count; i++)
        {
            distance = (target - utils::convertIdx1DTo3D(
                sampleData.data[i].index,
                sampleData.dimension
            )).magnitude();

            if(INVERSE)
            {
                interpolatedValue += C(i)/std::sqrtf(R*R + distance);
            }
            else
            {
                interpolatedValue += C(i)*std::sqrtf(R*R + distance);
            }
        }

        return { 
            originalIdx,
            interpolatedValue
        };
    }
}


const Eigen::MatrixXd precomputeGlobalC(const DataSource<ScatterPoint<float>>& sampleData, const float R)
{
    Eigen::MatrixXd M;
    Eigen::MatrixXd C;
    Eigen::MatrixXd F;

    M.resize(sampleData.count, sampleData.count);
    C.resize(sampleData.count);
    F.resize(sampleData.count);

    float distance;
    for(int i = 0; i < sampleData.count; i++)
    {
        F(i) = sampleData.data[i].funcValue;

        Point rowNeighbor = utils::convertIdx1DTo3D(
            sampleData.data[i].index,
            sampleData.dimension
        );

        for(int j = 0; j < sampleData.count; j++)
        {
            distance = (rowNeighbor - utils::convertIdx1DTo3D(
                sampleData.data[j].index,
                sampleData.dimension
            )).magnitude();

            M(i, j) = std::sqrtf(R*R + distance);
        }
    }

    C = M.householderQr().solve(F);

    return C;
}


}