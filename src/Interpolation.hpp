#pragma once


#include "Vector3.hpp"


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


}