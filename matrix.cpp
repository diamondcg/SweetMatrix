#include "matrix.h"
#include <math.h>

Mat4::Mat4(){
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

Mat4 Mat4::Rotate(double angle, Vec3 &axis){
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

Mat4 Mat4::Translate(Vec3 &t){
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

Mat4 Mat4::Scale(Vec3 &){
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

Mat4 Mat4::operator*(Vec4 &rh){

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

Mat4 Mat4::operator-(Mat4 &rh){
    Mat4 lh = *this;

    Mat4 m;
    for(int i=0; i<16; i++)
    {
        m.mat[i] = lh.mat[i]-rh.mat[i];
    }
    return m;
}

//Mat4& Mat4::operator-(double &rh){
//    Mat4 lh = *this;
//
//    Mat4 m;
//    m.mat = {
//         lh.mat[0]-rh,  lh.mat[1]-lh,  lh.mat[2]-rh,  lh.mat[3]-lh,
//         lh.mat[4]-rh,  lh.mat[5]-lh,  lh.mat[6]-rh,  lh.mat[7]-lh,
//         lh.mat[8]-rh,  lh.mat[9]-lh, lh.mat[10]-rh, lh.mat[11]-lh,
//        lh.mat[12]-rh, lh.mat[13]-lh, lh.mat[14]-rh, lh.mat[15]-lh
//    };
//
//    return m;
//}

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
