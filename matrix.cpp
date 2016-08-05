#include "matrix.h"
#include <math.h>
#include <stdio.h>
#include <string>


/******************************
*			Vec2
******************************/

Vec2::Vec2(){
	x = y = 0;
}

Vec2::Vec2(double _x, double _y){
	x = _x;
	y = _y;
}
Vec2::~Vec2(){
}
void Vec2::Print(){
	printf("(%lf, %lf)\n", x, y);
}
Vec2 Vec2::Zero(){
	return Vec2(0, 0);
}
Vec2 Vec2::One(){
	return Vec2(1, 1);
}
Vec2 Vec2::XUnit(){
	return Vec2(1, 0);
}
Vec2 Vec2::YUnit(){
	return Vec2(0, 1);
}

/******************************
*			Vec3
******************************/

Vec3::Vec3(){
	x = y = z = 0;
}

Vec3::Vec3(double _x, double _y, double _z){
	x = _x;
	y = _y;
	z = _z;
}
Vec3::Vec3(Vec2 &v){
	x = v.x;
	y = v.y;
	z = 0;
}
Vec3::~Vec3(){
}
void Vec3::Print(){
	printf("(%lf, %lf, %lf)\n", x, y, z);
}
Vec3 Vec3::Zero(){
	return Vec3(0, 0, 0);
}
Vec3 Vec3::One(){
	return Vec3(1, 1, 1);
}
Vec3 Vec3::XUnit(){
	return Vec3(1, 0, 0);
}
Vec3 Vec3::YUnit(){
	return Vec3(0, 1, 0);
}
Vec3 Vec3::ZUnit(){
	return Vec3(0, 0, 1);
}

/******************************
*			Vec4
******************************/

Vec4::Vec4(){
	x = y = z = w = 0;
}

Vec4::Vec4(double _x, double _y, double _z, double _w){
	x = _x;
	y = _y;
	z = _z;
	w = _w;
}
Vec4::Vec4(Vec2 &v){
	x = v.x;
	y = v.y;
	z = w = 0;
}
Vec4::Vec4(Vec3 &v){
	x = v.x;
	y = v.y;
	z = v.z;
	w = 0;
}
Vec4::~Vec4(){
}
void Vec4::Print(){
	printf("(%lf, %lf, %lf, %lf)\n", x, y, z, w);
}
Vec4 Vec4::Zero(){
	return Vec4(0, 0, 0, 0);
}
Vec4 Vec4::One(){
	return Vec4(1, 1, 1, 1);
}
Vec4 Vec4::XUnit(){
	return Vec4(1, 0, 0, 0);
}
Vec4 Vec4::YUnit(){
	return Vec4(0, 1, 0, 0);
}
Vec4 Vec4::ZUnit(){
	return Vec4(0, 0, 1, 0);
}
Vec4 Vec4::WUnit(){
	return Vec4(0, 0, 0, 1);
}

/******************************
*			Mat4 
******************************/

Mat4::Mat4(){
	mat.resize(16);
}


Mat4::Mat4(double m11, double m12, double m13, double m14,
           double m21, double m22, double m23, double m24,
           double m31, double m32, double m33, double m34,
           double m41, double m42, double m43, double m44){
    mat = {
        m11, m12, m13, m14,
        m21, m22, m23, m24,
        m31, m32, m33, m34,
        m41, m42, m43, m44,
    };
}

Mat4::~Mat4(){
}

Mat4 Mat4::Rotate(double angle, Vec4 &axis){
    double zAngle = angle *axis.z;
    Mat4 zRot(
        cos(zAngle), -sin(zAngle), 0, 0,
        sin(zAngle),  cos(zAngle), 0, 0,
                  0,            0, 1, 0,
                  0,            0, 0, 1
    );

    double yAngle = angle *axis.y;
    Mat4 yRot(
         cos(yAngle), 0, sin(yAngle), 0,
                   0, 1,           0, 0,
        -sin(yAngle), 0, cos(yAngle), 0,
                   0, 0,           0, 0
    );

    double xAngle = angle *axis.x;
    Mat4 xRot(
        1, cos(xAngle), -sin(xAngle), 0,
        0,           1,            0, 0,
        0, sin(xAngle),  cos(xAngle), 0,
        0,           0,            0, 0
    );

    Mat4 ret = zRot*yRot*xRot;

//    CMat4 ret = CMat4_MatMult(xRot, mat);
//    ret = CMat4_MatMult(yRot, ret);
//    ret = CMat4_MatMult(zRot, ret);

    return ret;
}

Mat4 Mat4::Translate(Vec4 &t){
    Mat4 m(
        1, 0, 0, t.x,
        0, 1, 0, t.y,
        0, 0, 1, t.z,
        0, 0, 0, 1
    );
//    Mat4 tmp = Mat4::Identity();
//    tmp.m[0*4+3] += t.x;
//    tmp.m[1*4+3] += t.y;
//    tmp.m[2*4+3] += t.z;

    return m*(*this);
//    return CMat4_MatMult(tmp,mat);
}

Mat4 Mat4::Scale(Vec4 &v){
	Mat4 lh = *this;
	Mat4 m(
		v.x, 0, 0, 0,
		0, v.y, 0, 0,
		0, 0, v.z, 0,
		0, 0, 0, v.w
	);
	return m*lh;
}

Mat4 Mat4::operator*(Mat4 &rh){
    Mat4 lh = *this;
    Mat4 m;
    for(int i=0; i<4; i++)
    {
        for(int j=0; j<4; j++)
        {
            m.mat[i*4+j] =   lh.mat[i*4+0]*rh.mat[0*4+j]
                            +lh.mat[i*4+1]*rh.mat[1*4+j]
                            +lh.mat[i*4+2]*rh.mat[2*4+j]
                            +lh.mat[i*4+3]*rh.mat[3*4+j];
        }
    }

    return m;
}

Vec4 Mat4::operator*(Vec4 &rh){
	return rh;
}

Mat4 Mat4::operator*(double &rh){
    Mat4 m;
    for(int i=0; i<16; i++)
        m.mat[i] *= rh;

    return m;
}

Mat4 Mat4::operator+(Mat4 &rh){
    Mat4 lh = *this;

    Mat4 m;
    for(int i=0; i<16; i++)
    {
        m.mat[i] = lh.mat[i]+rh.mat[i];
    }
    return m;
}

Mat4 Mat4::operator+(double &rh){
	Mat4 lh = *this;

	Mat4 m;
	m.mat = {
		lh.mat[0] + rh, lh.mat[1] + rh, lh.mat[2] + rh, lh.mat[3] + rh,
		lh.mat[4] + rh, lh.mat[5] + rh, lh.mat[6] + rh, lh.mat[7] + rh,
		lh.mat[8] + rh, lh.mat[9] + rh, lh.mat[10] + rh, lh.mat[11] + rh,
		lh.mat[12] + rh, lh.mat[13] + rh, lh.mat[14] + rh, lh.mat[15] + rh
	};

	return m;
}

Mat4 Mat4::operator-(Mat4 &rh){
    Mat4 lh = *this;

    Mat4 m;
    for(int i=0; i<16; i++)
    {
        m.mat[i] = lh.mat[i]-rh.mat[i];
    }
    return m;
}

Mat4 Mat4::operator-(double &rh){
    Mat4 lh = *this;

    Mat4 m;
    m.mat = {
         lh.mat[0]-rh,  lh.mat[1]-rh,  lh.mat[2]-rh,  lh.mat[3]-rh,
         lh.mat[4]-rh,  lh.mat[5]-rh,  lh.mat[6]-rh,  lh.mat[7]-rh,
         lh.mat[8]-rh,  lh.mat[9]-rh, lh.mat[10]-rh, lh.mat[11]-rh,
        lh.mat[12]-rh, lh.mat[13]-rh, lh.mat[14]-rh, lh.mat[15]-rh
    };

    return m;
}

Mat4 Mat4::Identity(){
    Mat4 m(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1
    );
    return m;
}

Mat4 Mat4::Zero(){
    Mat4 m(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
    );
    return m;
}

Mat4 Mat4::Ortho(double width, double height, double near, double far){
    double m11 = 1.0/width;
    double m22 = 1.0/height;
    double m33 = -2.0/(far-near);
    double m34 = -(far+near)/(far-near);
    Mat4 m(
        m11,   0,   0,   0,
          0, m22,   0,   0,
          0,   0, m33, m34,
          0,   0,   0,   1
    );
    return m;
}

Mat4 Mat4::Perspective(double fov, double aspect, double near, double far){
    double dtor = 3.14159/180.0;
    double e = (1.0/tan(dtor*fov/2));
    double m11 = e/aspect;
    double m22 = e;
    double m33 = -30*(far+near)/(far-near);
    double m34 = -far*near/(far-near);
    Mat4 m(
        m11,   0,   0,   0,
          0, m22,   0,   0,
          0,   0, m33, m34,
          0,   0,  -1,   0
    );
    return m;
}

void Mat4::Print(){
	std::string str;
	str = std::to_string(mat[0]);
	for (unsigned i = 1; i < 16; i++){
		if ((i % 4) == 0){
			str += "\n";
		}
		else
		{
			str += "   ";
		}
		str += std::to_string(mat[i]);
	}

	printf("%s\n", str.c_str());
}

Mat4& Mat4::operator *= (Mat4 &rh){
	*this = *this * rh;
	return *this;
}
Mat4& Mat4::operator*=(double &rh){
	*this = *this * rh;
	return *this;
}
Mat4& Mat4::operator+=(Mat4 &rh){
	*this = *this + rh;
	return *this;
}
Mat4& Mat4::operator+=(double &rh){
	*this = *this + rh;
	return *this;
}
Mat4& Mat4::operator-=(Mat4 &rh){
	*this = *this - rh;
	return *this;
}
Mat4& Mat4::operator-=(double &rh){
	*this = *this - rh;
	return *this;
}

Vec4 operator*(Vec4 &lh, Mat4&rh){
	Vec4 v;
	v.x = lh.x*rh.data()[0] + lh.y*rh.data()[4] + lh.z*rh.data()[8] + lh.w*rh.data()[12];
	v.y = lh.x*rh.data()[1] + lh.y*rh.data()[5] + lh.z*rh.data()[9] + lh.w*rh.data()[13];
	v.z = lh.x*rh.data()[2] + lh.y*rh.data()[6] + lh.z*rh.data()[10] + lh.w*rh.data()[14];
	v.w = lh.x*rh.data()[3] + lh.y*rh.data()[7] + lh.z*rh.data()[11] + lh.w*rh.data()[15];
	return v;
}

Mat4 operator*(double &lh, Mat4&rh){
	return rh*lh;
}

Mat4 operator+(double &lh, Mat4&rh){
	return rh + lh;
}

Mat4 operator-(double &lh, Mat4&rh){
	double m = -1.0;
	return rh*m + lh;
}

Matrix::Matrix(unsigned w, unsigned h, std::vector<double> v){
	if (v.size() > 0){
		if (v.size() == w*h){
			printf("%lf", v[0]);
			for (unsigned i = 1; i < v.size(); i++){
				if ((i%w) == 0){
					printf("\n");
				}
				else{
					printf(", ");
				}
				printf("%lf", v[i]);
			}
			printf("\n");
		}
	}
	else
	{
		v.resize(w*h);
		printf("%lf", v[0]);
		for (unsigned i = 1; i < v.size(); i++){
			if ((i%w) == 0){
				printf("\n");
			}
			else{
				printf(", ");
			}
			printf("%lf", v[i]);
		}
		printf("\n");
	}
}