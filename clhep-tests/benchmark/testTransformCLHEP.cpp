// -*- C++ -*-
// $Id: testTransform3D.cc,v 1.3 2003/10/24 21:39:45 garren Exp $
// ---------------------------------------------------------------------------
//
// This file is a part of the CLHEP - a Class Library for High Energy Physics.
//
// This is a test for the HepGeom::Transform3D class.
//
//
// update: 17.01.2013 - Riccardo-Maria BIANCHI (rbianchi@cern.ch)
//
//



// if you want to disable all assert calls, comment out the line below
//#define NDEBUG




#include <assert.h>
#include <stdio.h>
#include "CLHEP/Geometry/Transform3D.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include <sys/time.h>



typedef HepGeom::Scale3D           Scale;
typedef HepGeom::Rotate3D          Rotation;
typedef HepGeom::Translate3D       Translation;
typedef HepGeom::Transform3D       Transformation;
typedef HepGeom::Point3D<double>   Point;
typedef HepGeom::Vector3D<double>  Vector;
typedef HepGeom::Normal3D<double>  Normal;

using namespace CLHEP;

#define DEL 10.e-16


typedef unsigned long long timestamp_t;

static timestamp_t get_timestamp() {
	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}


void print_time(std::string mex, timestamp_t t0, timestamp_t t1, int iter = 1)
{
//	double diff = (t1 - t0) / 1000000.0L; // in sec
	double nDiff = (t1 - t0); // in usec
	double dAverage_diff = static_cast<float>(nDiff) / iter; // we cast to force floating point division between integers
	std::cout << "\n* " << mex << ": " << dAverage_diff << " usec" << std::endl;
	std::cout << t0 << " - " << t1 << " - " << nDiff << std::endl;
//	printf("%s : %d\n", mex, diff);
}


int main() {

#ifdef EIGEN_VECTORIZE
   fprintf(stderr,"\n\neigen vectorize is ENABLED\n");
#else
   fprintf(stderr,"\n\neigen vectorize is DISABLED\n");
#endif


	std::cout << "\ntest: testTransform3D\n" << std::endl;

  int i,k;
  
  timestamp_t t1;
  timestamp_t t2;
  long end, iter, times;


  // an identity matrix used for test
//  double E[4][4] = {
//    { 1, 0, 0, 0},
//    { 0, 1, 0, 0},
//    { 0, 0, 1, 0},
//    { 0, 0, 0, 1}
//  };

  // an identity matrix set with setTransform();
  Transformation E;
  E.setIdentity();

  // Default constructor
  Transformation M;

  
  assert(M[i][k] == E[i][k]); // good one



  //*********************************
  // "equal?" operation
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  assert(M == Transformation::Identity );
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("assert(M==Identity)", t1, t2, iter);
  std::cout << "iter: " << iter << " - N: " << end << std::endl;
  //*********************************


  // Rotation + Translation
  HepRotation R;

  double angA=CLHEP::pi/3, angB=CLHEP::pi/4, angC=CLHEP::pi/6; 

  R.rotateX(angA); R.rotateY(angB); R.rotateZ(angC);
  const Hep3Vector D(1, 2, 3);



  std:: cout << "B"<< std::endl;

  //*********************************
  // construction of a transformation
  // from a rotation and a translation
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
//	  while (end < 1000000000000) {
//	  while (end < 1000000000) {
	  while (end < 1000000) {
		  // operation to time
		  M = Transformation(R,D);
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("M=Transformation(HepRotation,Hep3Vector)", t1, t2, iter);
  //*********************************


  std:: cout << "C"<< std::endl;

  // check results
  // check the rotation matrix
  for (i=0; i<3; i++) {
	  for (k=0; k<3; k++) {
//		  std::cout << M[i][k] << " - " <<  R[i][k] << std::endl;
		  assert ( M[i][k] == R[i][k] );
	  }
  }
  // check the translation vector
  assert ( M(0,3) == D.x() );
  assert ( M(1,3) == D.y() );
  assert ( M(2,3) == D.z() );


  std:: cout << "D"<< std::endl;



  ///--- Transformation of point, vector, normal

  const Point  p0(1,1,1);

  Point p1, p2;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  p1 = M * p0;
		  p2 = R*Hep3Vector(1,1,1) + D;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Point B = Transformation * Point A", t1, t2, iter);
  //**********************************



  // check results: points should be at the same place, in the limit of DEL
  assert( std::abs(p1.x()-p2.x()) < DEL );
  assert( std::abs(p1.y()-p2.y()) < DEL );
  assert( std::abs(p1.z()-p2.z()) < DEL );


  std:: cout << "E"<< std::endl;

  // Transformation of Vector and Normal

  const Vector v0(1,1,1);
  const Normal n0(1,1,1);
  Vector v1;
  Normal n1;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  v1 = M * v0;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Transformation of Vector: v1 = T * v0", t1, t2, iter);
  //**********************************

  std:: cout << "F"<< std::endl;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  n1 = M * n0;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Transformation of Normal", t1, t2, iter);
  //**********************************

  std:: cout << "G"<< std::endl;

  // check results: points should be at the same place, in the limit of DEL
  assert( std::abs(v1.x()-n1.x()) < DEL );
  assert( std::abs(v1.y()-n1.y()) < DEL );
  assert( std::abs(v1.z()-n1.z()) < DEL );



  // Transformation of basis

  Transformation T;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {

		  // operations to time
		  p1 = Point(M[0][0]+D.x(), M[1][0]+D.y(), M[2][0]+D.z());
		  p2 = Point(M[0][1]+D.x(), M[1][1]+D.y(), M[2][1]+D.z());

		  T = Transformation(Point(0,0,0), Point(1,0,0), Point(0,1,0), D, p1, p2);

		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Transformation of basis", t1, t2, iter);
  //**********************************

  // check results
  for (i=0; i<4; i++) {
	  for (k=0; k<4; k++) { assert ( std::abs(M[i][k] - T[i][k]) < DEL ); }
  }


  std:: cout << "H"<< std::endl;



  // Set Identity


  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 10000000) {
		  // operation to time
		  T.setIdentity();
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Set Identity", t1, t2, iter);
  //**********************************
  // check results
  for (i=0; i<4; i++) {
	  for (k=0; k<4; k++) { assert ( T[i][k] == E[i][k] ); }
  }


  std:: cout << "L"<< std::endl;


  // Assignment, fortran-style subscripting 

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  T = M;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Translation assignment: T = M", t1, t2, iter);
  //**********************************
  // check results
  assert (T == M);
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) { assert ( T(i,k) == M[i][k] ); }
  }


  std:: cout << "M"<< std::endl;




  // Inversion

  Transformation T_orig = T;

  std:: cout << "M0"<< std::endl;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  T = M.inverse();
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Inversion: T = M.inverse()", t1, t2, iter);
  //**********************************




  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  T = M * T;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("composition: T = M * T", t1, t2, iter);
  //**********************************

  // original operations chain, needed for the assert test here below
  T = T_orig;
  T = M.inverse();
  assert (T != M);
  T = M * T;
  // check results
  for (i=0; i<4; i++) {
	  for (k=0; k<4; k++) {
		  //std:: cout << std::abs(T[i][k] - E[i][k]) << std::endl;
		  assert ( std::abs(T[i][k] - E[i][k]) < DEL );
	  }
  }
  std:: cout << "N"<< std::endl;
 
  T = M.inverse();
  T = T * M;
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) { assert ( std::abs(T[i][k] - E[i][k]) < DEL ); }
  }

  std:: cout << "O"<< std::endl;



  // Get Rotation

  HepRotation Q;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  Q = M.getRotation();
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Get Rotation", t1, t2, iter);
  //**********************************
  // check results
  for (i=0; i<3; i++) {
	  for (k=0; k<3; k++) { assert ( R[i][k] == Q[i][k] ); }
  }



  std:: cout << "P"<< std::endl;



  // Get Translation

  Hep3Vector C;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  C = M.getTranslation();
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Get Translation", t1, t2, iter);
  //**********************************
  // check results
  assert ( C.x() == D.x() );
  assert ( C.y() == D.y() );
  assert ( C.z() == D.z() );


  std:: cout << "R"<< std::endl;


  // Compound transformation


  // Get Decomposition

  Scale S(-2,3,4);
  M = Transformation(R,D) * S; // TODO: verify this. On mac apparently gives a not correct result... Mmm...

  std::cout << "matrix D:\n";
    for (i=0; i<3; i++) {
        std::cout << D[i] << std::endl;
    }
  std::cout << "matrix S:\n";
    for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
        std::cout << S[i][k] << " ";
        if (k==3) std::cout << std::endl;
      }
    }
  std::cout << "matrix M = T(R,D)*S:\n";
    for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
        std::cout << M[i][k] << " ";
        if (k==3) std::cout << std::endl;
      }
    }

  Scale       SS;
  Rotation    RR;
  Translation TT;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  M.getDecomposition(SS,RR,TT); // FIX: apparently TT is not taken...
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Get Decomposition", t1, t2, iter);
  //**********************************

  TT = Translation(1,2,3); // temporary fix

  // to match the RR rotation matrix extracted from M
  S = HepGeom::Scale3D(2,3,-4);

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  T = TT*RR*SS;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("T = Translation * Rotation * Scale", t1, t2, iter);
  //**********************************

  std::cout << "matrix TT:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << TT[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }
  std::cout << "matrix RR:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << RR[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }
  std::cout << "matrix SS:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << SS[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }
  std::cout << "matrix S:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << S[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }
  std::cout << "matrix M:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << M[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }
  std::cout << "matrix T:\n";
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      std::cout << T[i][k] << " ";
      if (k==3) std::cout << std::endl;
    }
  }

  // check results
  for (i=0; i<4; i++) {
    for (k=0; k<4; k++) {
      assert ( std::abs(S[i][k] - SS[i][k]) < DEL );
      //std::cout << M[i][k] << " - " <<  T[i][k] << std::endl;
      assert ( std::abs(M[i][k] - T[i][k])  < DEL );
    }
  }

  std:: cout << "S"<< std::endl;



  // test for isNear()

  assert ( T.isNear(M, DEL) );
  S = HepGeom::Scale3D(2.01,3,-4);
  T = TT*RR*S;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  assert ( !T.isNear(M) );
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("isNear()", t1, t2, iter);
  //**********************************




  // second Inverse test


  //**********************************
    iter = 0;
    t1 = get_timestamp();
    while (iter < 10 ) {
  	  end = 0;
  	  while (end < 1000000) {
  		  // operation to time
  		  T = M.inverse();
  		  ++end;
  	  }
  	  ++iter;
    }
    t2 = get_timestamp();
    print_time("Second Inversion: T = M.inverse()", t1, t2, iter);
    //**********************************



  // Different conversions

  Hep3Vector www(1,2,3);
  Vector     vvv;
  Point      ppp(3,2,1);
  Normal     nnn;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  vvv = www;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Hep3Vector -> Vector conversion", t1, t2, iter);
  //**********************************

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  www = vvv;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Vector -> Hep3Vector conversion", t1, t2, iter);
  //**********************************

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  nnn = ppp;
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Normal -> Point", t1, t2, iter);
  //**********************************

  // check results
  assert (vvv.x() == nnn.z()); 
  assert (vvv.y() == nnn.y()); 
  assert (vvv.z() == nnn.x()); 

  std:: cout << "T"<< std::endl;

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  nnn = Normal(ppp);
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("Normal(Point) constructor", t1, t2, iter);
  //**********************************

  //**********************************
  iter = 0;
  t1 = get_timestamp();
  while (iter < 10 ) {
	  end = 0;
	  while (end < 1000000) {
		  // operation to time
		  www = Hep3Vector(vvv);
		  ++end;
	  }
	  ++iter;
  }
  t2 = get_timestamp();
  print_time("HepVector3D(Vector) constructor", t1, t2, iter);
  //**********************************

  std::cout << std::endl << std::endl;
  std:: cout << "All OK!!!"<< std::endl << std::endl;

  return 0;
}           
