#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>

class Vec2{
protected:
    double x;
    double y;
public:
    Vec2();
    ~Vec2();
    static Vec2 Zero();
    static Vec2 One();
    static Vec2 XUnit();
    static Vec2 YUnit();
};

class Vec3{
public:
    double x;
    double y;
    double z;

    Vec3();
    ~Vec3();
    static Vec2 Zero();
    static Vec2 One();
    static Vec2 XUnit();
    static Vec2 YUnit();
    static Vec2 ZUnit();
};

class Vec4{
protected:
    double x;
    double y;
    double z;
    double w;
public:
    Vec4();
    Vec4(Vec2 &);
    Vec4(Vec3 &);
    ~Vec4();
    static Vec2 Zero();
    static Vec2 One();
    static Vec2 XUnit();
    static Vec2 YUnit();
    static Vec2 ZUnit();
    static Vec2 WUnit();
};

class Mat4{
protected:
    std::vector<double> mat;
public:
    Mat4();
    Mat4(double, double, double, double,
         double, double, double, double,
         double, double, double, double,
         double, double,  double, double);
    ~Mat4();
    Mat4 Rotate(double angle, Vec3 &);
    Mat4 Translate(Vec3 &);
    Mat4 Scale(Vec3 &);
    Mat4 operator*(Mat4 &rh);
    Mat4 operator*(Vec4 &rh);
    Mat4 operator*(double &rh);
    Mat4 operator+(Mat4 &rh);
//    Mat4& operator+(double &rh);
    Mat4 operator-(Mat4 &rh);
//    Mat4& operator-(double &rh);

    static Mat4 Identity();
    static Mat4 Zero();
    static Mat4 Ortho(double width, double height, double near, double far);
    static Mat4 Perspective(double fov, double aspect, double near, double far);
};

#endif // MATRIX_H_INCLUDED
