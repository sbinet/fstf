#include <immintrin.h>
#include <Eigen/Dense>
#include <CLHEP/Matrix/Matrix.h>
#include <CLHEP/Geometry/Transform3D.h>
#include <root/Math/SMatrix.h>

typedef ROOT::Math::SMatrix<double, 4> SMatrix44;
typedef ROOT::Math::SMatrix<double, 5, 5> SMatrix55;
typedef ROOT::Math::SMatrix<double, 5, 3> SMatrix53;
typedef ROOT::Math::SMatrix<double, 3, 5> SMatrix35;

using namespace Eigen;
using namespace CLHEP;

class FMatrixAVX{
  public:
    FMatrixAVX(){
      for(int i = 0; i < 16; i++){
	m_matrix[i] = i;
      }
    }

    void basicMul(FMatrixAVX &o, FMatrixAVX &tmp){
      for(int i = 0; i < 16; i+=4){
	for(int j = 0; j < 4; j++){
	  tmp.m_matrix[i+j] = m_matrix[i] * o.m_matrix[j] + \
		     m_matrix[i+1] * o.m_matrix[4 + j] + \
		     m_matrix[i+2] * o.m_matrix[8 + j] + \
		     m_matrix[i+3] * o.m_matrix[12 + j];
	}
      }
    }

    /* Unoptimized version 
    void optMul(FMatrixAVX &o, FMatrixAVX &tmp){
       __m256d a_line, b_line, r_line;

      for(int i = 0; i < 16; i+=4){
	a_line = _mm256_set1_pd(m_matrix[i]);
	b_line = _mm256_load_pd(&o.m_matrix[0]);
	r_line = _mm256_mul_pd(a_line, b_line);

	for(int j = 1; j < 4; j++){
	  a_line = _mm256_set1_pd(m_matrix[i+j]);
	  b_line = _mm256_load_pd(&o.m_matrix[j*4]);
	  r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line), r_line);
	}

	_mm256_store_pd(&tmp.m_matrix[i], r_line);
      }
    } */

    void optMul(FMatrixAVX &o, FMatrixAVX &tmp){
      __m256d a_line, b_line, r_line;

      __m256d b_line0 = _mm256_load_pd(&o.m_matrix[0]);
      __m256d b_line1 = _mm256_load_pd(&o.m_matrix[4]);
      __m256d b_line2 = _mm256_load_pd(&o.m_matrix[8]);
      __m256d b_line3 = _mm256_load_pd(&o.m_matrix[12]);

      a_line = _mm256_set1_pd(m_matrix[0]);
      r_line = _mm256_mul_pd(a_line, b_line0);

      a_line = _mm256_set1_pd(m_matrix[1]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line1), r_line);


      a_line = _mm256_set1_pd(m_matrix[2]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line2), r_line);

      a_line = _mm256_set1_pd(m_matrix[3]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line3), r_line);

      _mm256_store_pd(&tmp.m_matrix[0], r_line);


      a_line = _mm256_set1_pd(m_matrix[4]);
      r_line = _mm256_mul_pd(a_line, b_line0);

      a_line = _mm256_set1_pd(m_matrix[5]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line1), r_line);


      a_line = _mm256_set1_pd(m_matrix[6]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line2), r_line);

      a_line = _mm256_set1_pd(m_matrix[7]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line3), r_line);

      _mm256_store_pd(&tmp.m_matrix[4], r_line);


      a_line = _mm256_set1_pd(m_matrix[8]);
      r_line = _mm256_mul_pd(a_line, b_line0);

      a_line = _mm256_set1_pd(m_matrix[9]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line1), r_line);


      a_line = _mm256_set1_pd(m_matrix[10]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line2), r_line);

      a_line = _mm256_set1_pd(m_matrix[11]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line3), r_line);

      _mm256_store_pd(&tmp.m_matrix[8], r_line);


      a_line = _mm256_set1_pd(m_matrix[12]);
      r_line = _mm256_mul_pd(a_line, b_line0);

      a_line = _mm256_set1_pd(m_matrix[13]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line1), r_line);


      a_line = _mm256_set1_pd(m_matrix[14]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line2), r_line);

      a_line = _mm256_set1_pd(m_matrix[15]);
      r_line = _mm256_add_pd(_mm256_mul_pd(a_line, b_line3), r_line);

      _mm256_store_pd(&tmp.m_matrix[12], r_line);
    }

    __attribute__((aligned(32))) double m_matrix[16];
};

void timeBasicMatrixMultiplyAVX(FMatrixAVX &a, FMatrixAVX &b, FMatrixAVX &c);
void timeSIMDMatrixMultiplyAVX(FMatrixAVX &a, FMatrixAVX &b, FMatrixAVX &c);
void timeEigenMatrixMultiplyAVX(Matrix<double, 4, 4> &x, Matrix<double, 4, 4> &y, Matrix4d &z);
void timeCLHEPMatrixMultiply(HepMatrix &a, HepMatrix &b, HepMatrix &c);
void timeSMatrixMultiply(SMatrix44 &a, SMatrix44 &b, SMatrix44 &c);

void timeComplexCLHEP(CLHEP::HepMatrix &a, CLHEP::HepMatrix &b, CLHEP::HepMatrix &c, double alpha, double beta);
void timeComplexEigen(Matrix<double, 5, 3> &x, Matrix<double, 3, 5> &y, Matrix<double, 5, 5> &z, double alpha, double beta);
void timeComplexSMatrix(SMatrix53 &a, SMatrix35 &b, SMatrix55 &c, double alpha, double beta);
void timeMxNEigen(Matrix<double, 5, 3> &x, Matrix<double, 3, 5> &y, Matrix<double, 5, 5> &z);
void timeMxNSMatrix(SMatrix53 &a, SMatrix35 &b, SMatrix55 &c);
