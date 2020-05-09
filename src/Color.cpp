#include <iostream>
#include <cmath>

#include "Color.hpp"


Color::Color(): _data{1.0f, 1.0f, 1.0f}
{
}

Color::Color(const float &f): _data{f, f, f}
{
}

Color::Color(const float &r, const float &g, const float &b): _data{r, g, b}
{
}

Color::Color(const Color &c): _data{c.r(), c.g(), c.b()}
{
}

Color::~Color()
{
}


Color Color::operator - () const
{
    return Color(-_data[0], -_data[1], -_data[2]);
}

float Color::operator [] (const int &idx) const
{
    return _data[idx];
}

float& Color::operator [] (const int &idx)
{
    return _data[idx];
}

Color& Color::operator = (const Color &c)
{
    _data[0] = c.r();
    _data[1] = c.g();
    _data[2] = c.b();

    return *this;
}

Color& Color::operator += (const Color &c)
{
    _data[0] += c.r();
    _data[1] += c.g();
    _data[2] += c.b();

    return *this;
}

Color& Color::operator += (const float &f)
{
    _data[0] += f;
    _data[1] += f;
    _data[2] += f;

    return *this;
}

Color& Color::operator -= (const Color &c)
{
    _data[0] -= c.r();
    _data[1] -= c.g();
    _data[2] -= c.b();

    return *this;
}

Color& Color::operator -= (const float &f)
{
    _data[0] -= f;
    _data[1] -= f;
    _data[2] -= f;

    return *this;
}

Color& Color::operator *= (const Color &c)
{
    _data[0] *= c.r();
    _data[1] *= c.g();
    _data[2] *= c.b();

    return *this;
}

Color& Color::operator *= (const float &f)
{
    _data[0] *= f;
    _data[1] *= f;
    _data[2] *= f;

    return *this;
}

Color& Color::operator *= (const Vector3 &v)
{
    _data[0] *= v.x();
    _data[1] *= v.y();
    _data[2] *= v.z();

    return *this;
}

Color& Color::operator /= (const Color &c)
{
    _data[0] /= c.r();
    _data[1] /= c.g();
    _data[2] /= c.b();

    return *this;
}

Color& Color::operator /= (const float &f)
{
    _data[0] /= f;
    _data[1] /= f;
    _data[2] /= f;

    return *this;
}


float Color::r() const
{
    return _data[0];
}

float Color::g() const
{
    return _data[1];
}

float Color::b() const
{
    return _data[2];
}

float& Color::r_()
{
    return _data[0];
}

float& Color::y_()
{
    return _data[1];
}

float& Color::z_()
{
    return _data[2];
}

void Color::setR(const float &f)
{
    _data[0] = f;
}

void Color::setG(const float &f)
{
    _data[1] = f;
}

void Color::setB(const float &f)
{
    _data[2] = f;
}

void Color::set(const float &r, const float &g, const float &b)
{
    _data[0] = r;
    _data[1] = g;
    _data[2] = b;
}

void Color::set(const float &f)
{
    _data[0] = f;
    _data[1] = f;
    _data[2] = f;
}

Color Color::clamp()
{
    float r = (_data[0] < 0.0f)? 0.0f : _data[0];
    r = (r > 1.0f)? 1.0f : r;

    float g = (_data[1] < 0.0f)? 0.0f : _data[1];
    g = (g > 1.0f)? 1.0f : g;

    float b = (_data[2] < 0.0f)? 0.0f : _data[2];
    b = (b > 1.0f)? 1.0f : b;

    return { r, g, b };
}

Color Color::convert255() const
{
    return Color(int(255.99*_data[0]), int(255.99*_data[1]), int(255.99*_data[2]));
}

Color Color::gammaCorrection(float gamma) const 
{
    float invGamma = 1/gamma;
    return Color(std::pow(_data[0], invGamma), std::pow(_data[1], invGamma), std::pow(_data[2], invGamma));
}

Color Color::gammaCorrection2() const
{
    return Color(std::sqrt(_data[0]), std::sqrt(_data[1]), std::sqrt(_data[2]));
}


Color operator + (const Color &c1, const Color &c2)
{
    return Color(c1.r() + c2.r(), c1.g() + c2.g(), c1.b() + c2.b());
}

Color operator + (const Color &c, const float &f)
{
    return Color(c.r() + f, c.g() + f, c.b() + f);
}

Color operator + (const float &f, const Color &c)
{
    return Color(c.r() + f, c.g() + f, c.b() + f);
}

Color operator - (const Color &c1, const Color &c2)
{
    return Color(c1.r() - c2.r(), c1.g() - c2.g(), c1.b() - c2.b());
}

Color operator - (const Color &c, const float &f)
{
    return Color(c.r() - f, c.g() - f, c.b() - f);
}

Color operator - (const float &f, const Color &c)
{
    return Color(f - c.r(), f - c.g(), f - c.b());
}

Color operator * (const Color &c1, const Color &c2)
{
    return Color(c1.r()*c2.r(), c1.g()*c2.g(), c1.b()*c2.b());
}

Color operator * (const Color &c, const float &f)
{
    return Color(c.r()*f, c.g()*f, c.b()*f);
}

Color operator * (const float &f, const Color &c)
{
    return Color(c.r()*f, c.g()*f, c.b()*f);
}

Color operator * (const Vector3 &v, const Color &c)
{
    return Color(v.x()*c.r(), v.y()*c.g(), v.z()*c.b());
}

Color operator * (const Color &c, const Vector3 &v)
{
    return Color(v.x()*c.r(), v.y()*c.g(), v.z()*c.b());
}

Color operator / (const Color &c1, const Color &c2) 
{
    return Color(c1.r()/c2.r(), c1.g()/c2.g(), c1.b()/c2.b());
}

Color operator / (const Color &c, const float &f) 
{
    return Color(c.r()/f, c.g()/f, c.b()/f);
}

Color operator / (const float &f, const Color &c)
{
    return Color(f/c.r(), f/c.g(), f/c.b());
}

std::istream &operator >> (std::istream &s, Color &c)
{
    s >> c.r_() >> c.y_() >> c.z_();

    return s;
}

std::ostream &operator << (std::ostream &s, const Color &c)
{
    s << c.r() << " " << c.g() << " " << c.b();

    return s;
}
