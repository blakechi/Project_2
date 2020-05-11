#pragma once

template<typename T>
struct Vector2
{
    union
    {
        struct 
        {
            T x, y;
        };
        
        struct
        {
            T min, max;
        };
    };
};

template<typename T>
using Point2 = Vector2<T>;