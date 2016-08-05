#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <vector>

class Vec2{
public:
	double x;
	double y;
	Vec2();
	Vec2(double _x, double _y);
	~Vec2();
	void Print();
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
	Vec3(double _x, double _y, double _z);
	Vec3(Vec2 &);
	~Vec3();
	void Print();
    static Vec3 Zero();
    static Vec3 One();
    static Vec3 XUnit();
    static Vec3 YUnit();
    static Vec3 ZUnit();
};

class Vec4{
public:
	double x;
	double y;
	double z;
	double w;
    Vec4();
	Vec4(double _x, double _y, double _z, double _w);
    Vec4(Vec2 &);
    Vec4(Vec3 &);
    ~Vec4();
	void Print();
    static Vec4 Zero();
    static Vec4 One();
    static Vec4 XUnit();
    static Vec4 YUnit();
    static Vec4 ZUnit();
    static Vec4 WUnit();
};

class Matrix{
public:
	Matrix(unsigned w, unsigned h, std::vector<double> v = std::vector<double>());
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
    Mat4 Rotate(double angle, Vec4 &);
    Mat4 Translate(Vec4 &);
    Mat4 Scale(Vec4 &);
    Mat4 operator*(Mat4 &rh);
    Vec4 operator*(Vec4 &rh);
    Mat4 operator*(double &rh);
    Mat4 operator+(Mat4 &rh);
    Mat4 operator+(double &rh);
    Mat4 operator-(Mat4 &rh);
	Mat4 operator-(double &rh);
	Mat4& operator*=(Mat4 &rh);
	Mat4& operator*=(double &rh);
	Mat4& operator+=(Mat4 &rh);
	Mat4& operator+=(double &rh);
	Mat4& operator-=(Mat4 &rh);
	Mat4& operator-=(double &rh);
	const double* data() { return mat.data(); }
	void Print();

    static Mat4 Identity();
    static Mat4 Zero();
    static Mat4 Ortho(double width, double height, double near, double far);
    static Mat4 Perspective(double fov, double aspect, double near, double far);
};

Vec4 operator*(Vec4 &lh, Mat4&rh);
Mat4 operator*(double &lh, Mat4&rh);
Mat4 operator+(double &lh, Mat4&rh);
Mat4 operator-(double &lh, Mat4&rh);

#endif // MATRIX_H_INCLUDED
