#include "testAVX.h"

void timeBasicMatrixMultiplyAVX(FMatrixAVX &a, FMatrixAVX &b, FMatrixAVX &c){
  a.basicMul(b, c);
}

void timeSIMDMatrixMultiplyAVX(FMatrixAVX &a, FMatrixAVX &b, FMatrixAVX &c){
  a.optMul(b, c);
}

void timeEigenMatrixMultiplyAVX(Matrix<double, 4, 4> &x, Matrix<double, 4, 4> &y, Matrix4d &z){
  z.noalias() = x * y;
}

void timeCLHEPMatrixMultiply(HepMatrix &a, HepMatrix &b, HepMatrix &c){
  c = a * b;
}

void timeSMatrixMultiply(SMatrix44 &a, SMatrix44 &b, SMatrix44 &c){
  c = a * b;
}

void timeComplexCLHEP(HepMatrix &a, HepMatrix &b, HepMatrix &c, double alpha, double beta){
  c = alpha * (a * b) + beta * c;
}

void timeComplexEigen(Matrix<double, 5, 3> &a, Matrix<double, 3, 5> &b, Matrix<double, 5, 5> &c, double alpha, double beta){
  c = alpha * (a * b) + beta * c;
}

void timeComplexSMatrix(SMatrix53 &a, SMatrix35 &b, SMatrix55 &c, double alpha, double beta){
  c = alpha * (a * b) + beta * c;
}

void timeMxNEigen(Matrix<double, 5, 3> &a, Matrix<double, 3, 5> &b, Matrix<double, 5, 5> &c){
  c.noalias() = a * b;
}

void timeMxNSMatrix(SMatrix53 &a, SMatrix35 &b, SMatrix55 &c){
  c = a * b;
}
