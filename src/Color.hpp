#pragma once


#include "Vector3.hpp"


class Color
{
private:
    float _data[3];

public:
    Color();
    Color(const float &f);
    Color(const float &r, const float &g, const float &b);
    Color(const Color &c);
    ~Color();

    Color operator - () const;
    float operator [] (const int &idx) const;
    float& operator [] (const int &idx);
    Color& operator = (const Color &c);
    Color& operator += (const Color &c);
    Color& operator += (const float &f);
    Color& operator -= (const Color &c);
    Color& operator -= (const float &f);
    Color& operator *= (const Color &c);
    Color& operator *= (const float &f);
    Color& operator *= (const Vector3 &v);
    Color& operator /= (const Color &c);
    Color& operator /= (const float &f);

    float r() const;
    float g() const;
    float b() const;
    float& r_();
    float& y_();
    float& z_();
    void setR(const float &f);
    void setG(const float &f);
    void setB(const float &f);
    void set(const float &r, const float &g, const float &b);
    void set(const float &f);
    Color clamp();
    Color convert255() const;
    Color gammaCorrection(float gamma) const;
    Color gammaCorrection2() const;
};

Color operator + (const Color &c1, const Color &c2);
Color operator + (const Color &c, const float &f);
Color operator + (const float &f, const Color &c);
Color operator - (const Color &c1, const Color &c2);
Color operator - (const Color &c, const float &f);
Color operator - (const float &f, const Color &c);
Color operator * (const Color &c1, const Color &c2);
Color operator * (const Color &c, const float &f);
Color operator * (const float &f, const Color &c);
Color operator * (const Vector3 &v, const Color &c);
Color operator * (const Color &c, const Vector3 &v);
Color operator / (const Color &c1, const Color &c2);
Color operator / (const Color &c, const float &f);
Color operator / (const float &f, const Color &c);
std::istream &operator >> (std::istream &s, Color &c);
std::ostream &operator << (std::ostream &s, const Color &c);

typedef Color Image;
