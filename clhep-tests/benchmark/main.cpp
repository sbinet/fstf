#include <iostream>
#include <sys/time.h>
#include <unistd.h>
#include <omp.h>
#include <vector>
#include <xmmintrin.h>

#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Geometry/Transform3D.h>
#include <root/Math/SMatrix.h>
#include <mkl.h>
#include <ipp.h>

#include "testAVX.h"

using namespace std;

double start, end;

void ippinfo(void){
  const IppLibraryVersion* lib = ippsGetLibVersion();
  cout << lib->Name << lib->Version << lib->Version << lib->major << lib->majorBuild << lib->build << endl;
}

void start_clock(){
  start = omp_get_wtime();
}

double stop_clock(){
  return omp_get_wtime() - start;
}

#define N 1000000

void initializeHEPMatrix(CLHEP::HepMatrix &hepM, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      hepM[i][j] = i * n + j;
    }
  }
}

void printMatrixAVX(FMatrixAVX &a){
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      cout << a.m_matrix[i*4 + j] << " ";
    }
    cout << endl;
  }
}

void printHEPMatrix(CLHEP::HepMatrix &a, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      cout << a[i][j] << " ";
    }
    cout << endl;
  }
}

void printdgemm(double *a, int m, int n){
  for(int i = 0; i < m; i++){
    for(int j = 0; j < n; j++){
      cout << a[i*m + j] << " ";
    }
    cout << endl;
  }
}

void printSMatrix(SMatrix44 &a){
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      cout << a(i,j) << " ";
    }
    cout << endl;
  }
}

void printSMatrix(SMatrix55 &a){
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      cout << a(i,j) << " ";
    }

    cout << endl;
  }
}

double CLHEP_multiply(){
  CLHEP::HepMatrix a(4, 4), b(4, 4), c(4, 4);
  double timing;

  initializeHEPMatrix(a, 4, 4);
  initializeHEPMatrix(b, 4, 4);
  initializeHEPMatrix(c, 4, 4);

  start_clock();
  for(int i = 0; i < N; i++){
    timeCLHEPMatrixMultiply(a, b, c);
  }
  timing = stop_clock();
  cout << "CLHEP: " << timing << endl;

  printHEPMatrix(c, 4, 4);
  return timing;
}

void SIMD_multiply(double clhep){
  FMatrixAVX a, b, c;
  double timing;

  start_clock();
  for(int i = 0; i < N; i++){
    timeSIMDMatrixMultiplyAVX(a, b, c);
  }
  timing = stop_clock();
  cout << "AVX SIMD: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  printMatrixAVX(c);
}

void Basic_multiply(double clhep){
  FMatrixAVX a, b, c;
  double timing;

  start_clock();
  for(int i = 0; i < N; i++){
    timeBasicMatrixMultiplyAVX(a, b, c);
  }
  timing = stop_clock();
  cout << "Basic:" << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;


  printMatrixAVX(c);
}

void Eigen_multiply(double clhep){
  Matrix4d a, b, c;
  double timing;

  a << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;
  b << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15;

  start_clock();
  for(int i = 0; i < N; i++){
    timeEigenMatrixMultiplyAVX(a, b, c);
  }
  timing = stop_clock();
  cout << "Eigen: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  cout << c << endl;
}

void MKL_multiply(double clhep){
  double alpha = 1.0, beta = 0., timing;
  int n = 4;
  double a[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  double b[16] ={0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  double c[16];

  start_clock();
  for(int i = 0; i < N; i++){
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,alpha,b,n,a,n,beta,c,n);
  }
  timing = stop_clock();
  cout << "MKL: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  printdgemm(c, 4, 4);
}

void SMatrix_multiply(double clhep){
  double timing;
  double ar[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  SMatrix44 a(ar, 16), b(ar, 16), c(ar, 16);

  start_clock();
  for(int i = 0; i < N; i++){
    timeSMatrixMultiply(a, b, c);
  }
  timing = stop_clock();
  cout << "SMatrix: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  printSMatrix(c);
}

void IPP_multiply(double clhep){
  double timing;
  Ipp64f a[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Ipp64f b[16] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  Ipp64f c[16];

  int stride1 = sizeof(Ipp64f)*4;
  int stride2 = sizeof(Ipp64f);

  start_clock();
  for(int i = 0; i < N; i++)
    ippmMul_mm_64f(a, stride1, stride2, 4, 4, b, stride1, stride2, 4, 4, c, stride1, stride2);
  timing = stop_clock();
  cout << "IPP: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;
  ippinfo();

  printdgemm(c, 4, 4);
}

double CLHEP_complex(){
  HepMatrix a(5, 3), b(3, 5), c(5, 5);
  double timing;

  initializeHEPMatrix(a, 5, 3);
  initializeHEPMatrix(b, 3, 5);
  initializeHEPMatrix(c, 5, 5);

  start_clock();
  for(int i = 0; i < N; i++){
    timeComplexCLHEP(a, b, c, 1, 1);
  }
  timing = stop_clock();

  cout << "Complex CLHEP: " << timing << endl;
  printHEPMatrix(c, 5, 5);
  return timing;
}

void MKL_complex(double clhep, double alpha, double beta){
  double a[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  double b[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  double c[25] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  double timing;

  start_clock();
  for(int i = 0; i < N; i++){
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,5,5,3,alpha,a,3,b,5,beta,c,5);
  }
  timing = stop_clock();
  cout << "Complex/MxN MKL: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  printdgemm(c, 5, 5);
}

void Eigen_complex(double clhep){
  Matrix<double, 5, 3> a;
  Matrix<double, 3, 5> b;
  Matrix<double, 5, 5> c;
  double timing;

  a << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
  b << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
  c << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24;

  start_clock();
  for(int i = 0; i < N; i++){
    timeComplexEigen(a, b, c, 1, 1);
  }
  timing = stop_clock();

  cout << "Complex Eigen: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;
  cout << c << endl;
}

void SMatrix_complex(double clhep){
  double ar[25] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
  double timing;
  SMatrix53 a(ar, 15);
  SMatrix35 b(ar, 15);
  SMatrix55 c(ar, 25);

  for(int i = 0; i < N; i++){
    timeComplexSMatrix(a, b, c, 1, 1);
  }
  timing = stop_clock();

  cout << "Complex SMatrix: " <<  timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;  
  printSMatrix(c);
}

double CLHEP_MxN(){
  HepMatrix a(5, 3), b(3, 5), c(5, 5);
  double timing;

  initializeHEPMatrix(a, 5, 3);
  initializeHEPMatrix(b, 3, 5);
  initializeHEPMatrix(c, 5, 5);

  start_clock();
  for(int i = 0; i < N; i++){
    timeCLHEPMatrixMultiply(a, b, c);
  }
  timing = stop_clock();

  cout << "MxN CLHEP: " <<  endl;
  printHEPMatrix(c, 5, 5);
  return timing;
}

double Eigen_MxN(double clhep){
  Matrix<double, 5, 3> a;
  Matrix<double, 3, 5> b;
  Matrix<double, 5, 5> c;
  double timing;

  a << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;
  b << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14;

  start_clock();
  for(int i = 0; i < N; i++){
    timeMxNEigen(a, b, c);
  }
  timing = stop_clock();

  cout << "MxN Eigen: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;
  cout << c << endl;
}

double SMatrix_MxN(double clhep){
  double ar[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  double timing;
  SMatrix53 a(ar, 15);
  SMatrix35 b(ar, 15);
  SMatrix55 c(ar, 15);

  for(int i = 0; i < N; i++){
    timeMxNSMatrix(a, b, c);
  }
  timing = stop_clock();

  cout << "MxN SMatrix: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;
  printSMatrix(c);
}

void IPP_MxN(double clhep){
  double timing;
  Ipp64f a[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  Ipp64f b[15] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
  Ipp64f c[25];

  int stride1 = sizeof(Ipp64f)*4;
  int stride2 = sizeof(Ipp64f);

  start_clock();
  for(int i = 0; i < N; i++)
    ippmMul_mm_64f(a, stride1, stride2, 5, 3, b, stride1, stride2, 3, 5, c, stride1, stride2);
  timing = stop_clock();
  cout << "MxN IPP: " << timing << "\nSpeedup vs CLHEP: " << clhep/timing << endl;

  printdgemm(c, 5, 5);
}

int main(){
  double clhep = CLHEP_multiply();
  Basic_multiply(clhep);
  SIMD_multiply(clhep);
  MKL_multiply(clhep);
  Eigen_multiply(clhep);
  SMatrix_multiply(clhep);
  IPP_multiply(clhep);

  cout << "\n\n";

  clhep = CLHEP_MxN();
  MKL_complex(clhep, 1, 0);
  Eigen_MxN(clhep);
  SMatrix_MxN(clhep);
  IPP_MxN(clhep);

  cout << "\n\n";

  clhep = CLHEP_complex();
  MKL_complex(clhep, 1, 1);
  Eigen_complex(clhep);
  SMatrix_complex(clhep);
}
