#pragma once

template<typename T>
struct DataSource
{
    T* data;
    int* dimension; // [x, y, z]
    int count;
    int maxValue;
    int minValue;
};