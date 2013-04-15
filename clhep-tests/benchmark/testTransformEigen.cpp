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
#include <iostream>
#include <sys/time.h>

#include <Eigen/Geometry>

#include "CLHEP/Units/PhysicalConstants.h"
#define DEL 10.e-15

typedef Eigen::Affine3d Transformation; 
typedef Eigen::Vector3d Point;
typedef Eigen::Vector3d Vector;
typedef Eigen::Translation3d Translation;

typedef unsigned long long timestamp_t;
static timestamp_t get_timestamp() {
	struct timeval now;
	gettimeofday(&now, NULL);
	return now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}

void print_time(std::string mex, timestamp_t t0, timestamp_t t1, int iter = 1)
{
	double nDiff = (t1 - t0); // in usec
	double dAverage_diff = static_cast<float>(nDiff) / iter;
	std::cout << "\n* " << mex << ": " << dAverage_diff << " usec" << std::endl;
}

int main() {
	int i,k;
	long end, iter;
	timestamp_t t1;
	timestamp_t t2;

	double E[4][4] = {
		{ 1, 0, 0, 0},
		{ 0, 1, 0, 0},
		{ 0, 0, 1, 0},
		{ 0, 0, 0, 1}
	};

	std::cout << "\ntest: testTransformationEigen\n" << std::endl;

	Transformation M = Transformation::Identity();

	// Rotation + Translation
	double angA=CLHEP::pi/3, angB=CLHEP::pi/4, angC=CLHEP::pi/6; 

	Transformation R = Eigen::AngleAxisd(angC, Vector::UnitZ()) * 
	    		   Eigen::AngleAxisd(angB, Vector::UnitY()) * 
	    		   Eigen::AngleAxisd(angA, Vector::UnitX()) *
			   Transformation::Identity();

	Translation D(1, 2, 3);

	//*********************************
	// construction of a transformation
	// from a rotation and a translation
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			M = D * R;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("M=Transformation(R,D)", t1, t2, iter);
	//*********************************

	std:: cout << "C"<< std::endl;

	std::cout << "Rotated matrix - R: \n" << R.matrix() << std::endl;
	std::cout << "default Transformation - M: \n" << M.matrix() << std::endl;

	// check the rotation matrix
	for (i=0; i<3; i++) {
		for (k=0; k<3; k++) {
			assert ( M(i,k) == R(i,k) );
		}
	}

	// check the translation vector
	assert ( M(0,3) == D.x() );
	assert ( M(1,3) == D.y() );
	assert ( M(2,3) == D.z() );

	std::cout << "D"<< std::endl;

	// Transformation of point, vector, normal
	Point p0(1,1,1), p1, p2, d = D.vector();

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			p1 = M * p0;
			p2 = R.linear() * Vector(1,1,1) + d;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Point B = Transformation * Point A", t1, t2, iter);
	//**********************************

	// check results: points should be at the same place, in the limit of DEL
	double xp_diff = p1(0)-p2(0);
	double yp_diff = p1(1)-p2(1);
	double zp_diff = p1(2)-p2(2);

	assert( std::abs(xp_diff) < DEL );
	assert( std::abs(yp_diff) < DEL );
	assert( std::abs(zp_diff) < DEL );

	std:: cout << "E"<< std::endl;

	// Transformation of Vector and Normal
	Vector v0(1,1,1), n0(1,1,1), v1, n1;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			v1 = M * v0;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Transformation of Vector", t1, t2, iter);
	//**********************************

	std:: cout << "F"<< std::endl;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
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
	double x_diff = v1(0)-n1(0);
	double y_diff = v1(1)-n1(1);
	double z_diff = v1(2)-n1(2);
	assert( std::abs(x_diff) < DEL );
	assert( std::abs(y_diff) < DEL );
	assert( std::abs(z_diff) < DEL );

	// Transformation of basis
	Transformation T, T2, T3 = Transformation::Identity();
	Vector x1,y1,z1, x2,y2,z2;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000 ) {
			Point to0 = D.vector();
			Point to1 = Point( M(0,0)+D.x(), M(1,0)+D.y(), M(2,0)+ D.z() );
			Point to2 = Point( M(0,1)+D.x(), M(1,1)+D.y(), M(2,1)+ D.z() );

			Point fr0 = Point(0,0,0);
			Point fr1 = Point(1,0,0);
			Point fr2 = Point(0,1,0);

			// Axes of the coordinate system "fr"
			x1 = (fr1 - fr0).normalized();
			y1 = (fr2 - fr0).normalized();

			// Axes of the coordinate system "to"
			x2 = (to1 - to0).normalized();
			y2 = (to2 - to0).normalized();

			// transformation matrix
			/*
			 * From:
			 * - http://gamedev.stackexchange.com/questions/26084/transform-between-two-3d-cartesian-coordinate-systems
			 * - http://stackoverflow.com/questions/15252919/how-to-find-the-transformation-matrix-of-a-change-of-basis-with-eigen
			 */

			T2.linear() << x1, y1, x1.cross(y1); // transform from CS1 to CS2 (in our simple case CS1 == CS2, so T2=Identity
			T3.linear() << x2, y2, x2.cross(y2); // transform from CS1 to CS3

			T.linear() = T3.linear() * T2.linear().inverse(); // T = transform to CS2 to CS3 (in our simple case: T = T3)
			T.translation() = D.vector();

			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Transformation of basis", t1, t2, iter);
	//**********************************

	// check results
	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) {
			double diff = M(i,k) - T(i,k);
			assert ( std::abs(diff) < DEL );
		}
	}

	std:: cout << "H"<< std::endl;

	// Set Identity

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 10000000) {
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
		for (k=0; k<4; k++) { assert ( T(i,k) == E[i][k] ); }
	}

	std:: cout << "L"<< std::endl;

	// Assignment, fortran-style subscripting 

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			T = M;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Translation assignment", t1, t2, iter);

	//**********************************
	// check results
	//  assert (T == M);
	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) { assert ( T(i,k) == M(i,k) ); }
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
			T = M.inverse();
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Inversion", t1, t2, iter);
	//**********************************

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			// operation to time
			T = M.inverse(Eigen::Affine);
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Inversion(Affine)", t1, t2, iter);
	//**********************************

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			// operation to time
			T = M.inverse(Eigen::Isometry);
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Inversion(Isometry)", t1, t2, iter);
	//**********************************

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			T = M * T;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("composition T = Transformation * T", t1, t2, iter);
	//**********************************

	// original operations chain, needed for the assert test here below
	T = T_orig;
	T = M.inverse();
	bool testOk = false;

	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) {
			//std::cout << "i-k:" << i << k << " - " << T(i,k) << " " << M(i,k) << std::endl;
			if (T(i,k) != M(i,k))  {
				testOk = true;
			}
		}
	}
	if (testOk==false) assert(false);

	T = M * T;
	// check results
	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) {
			double diff = T(i,k) - E[i][k];
			assert ( std::abs(diff) < DEL );
		}
	}
	std:: cout << "N"<< std::endl;

	T = M.inverse();
	T = T * M;

	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) {
			double diff = T(i,k) - E[i][k];
			assert ( std::abs(diff) < DEL );
		}
	}

	std:: cout << "O"<< std::endl;

	// Get Rotation
	Transformation Q;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			Q.linear() = M.linear();
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Get Rotation", t1, t2, iter);
	//**********************************
	
	// check results
	for (i=0; i<3; i++) {
		for (k=0; k<3; k++) {
			double a1 = R(i,k);
			double a2 = Q(i,k);
			assert ( a1 - a2 < DEL );
		}
	}

	std:: cout << "P"<< std::endl;

	// Get Translation
	Vector C;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			C = M.translation();
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Get Translation", t1, t2, iter);
	//**********************************
	// check results
	assert ( C(0) == D.x() );
	assert ( C(1) == D.y() );
	assert ( C(2) == D.z() );

	std:: cout << "R"<< std::endl;

	// Compound transformation
	// Get Decomposition
	Transformation S = Transformation::Identity();
	Vector VV(-2,3,4);
	S.scale(VV);

	Transformation M2;
	M2.linear() = R.linear();
	M2.translation() = D.vector();

	M = M2 * S;
	Transformation SS = Transformation::Identity();
	Transformation RR = Transformation::Identity();
	Translation TT = Translation::Identity();

	Eigen::Matrix3d mat_rotation, mat_scaling;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			// decomposing the linear part of M in a product of a rotation matrix and a scaling matrix
			M.computeRotationScaling(&mat_rotation, &mat_scaling);
			RR.linear() = mat_rotation;
			SS.linear() = mat_scaling;

			// get the translation part of M
			TT.translation() = M.translation();
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Get Decomposition", t1, t2, iter);
	//**********************************

	//check
	Eigen::Matrix3d md = mat_rotation * mat_scaling;
	assert( M.linear().isApprox(md, DEL) );

	// we set S to (2,3,-4) to match the scale part we got from M here above
	S = Transformation::Identity();
	Vector VVV(2,3,-4);
	S.scale(VVV);

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			// operation to time
			T = TT * RR * SS;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("T = Translation * Rotation * Scale", t1, t2, iter);
	//**********************************

	std::cout << "TT: \n" << TT.vector() << std::endl;
	std::cout << "RR: \n" << RR.matrix() << std::endl;
	std::cout << "SS: \n" << SS.matrix() << std::endl;
	std::cout << "S: \n" << S.matrix() << std::endl;
	std::cout << "M: \n" << M.matrix() << std::endl;
	std::cout << "T: \n" << T.matrix() << std::endl;

	// check results
	for (i=0; i<4; i++) {
		for (k=0; k<4; k++) {
			double diff_a = S(i,k) - SS(i,k);
			double diff_b = M(i,k) - T(i,k);
			//std::cout << "a-b: " << diff_a << " " << diff_b << std::endl;
			assert ( std::abs(diff_a) < DEL );
			assert ( std::abs(diff_b)  < DEL );
		}
	}

	std:: cout << "S"<< std::endl;

	assert ( T.isApprox(M, DEL) );

	Transformation S2 = Transformation::Identity();
	Vector VV2(2.01 ,3, -4);

	S2.scale(VV2);
	T = TT * RR * S2;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			assert ( !T.isApprox(M) );
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
			T = M.inverse(Eigen::Affine);
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Second Inversion: T = M.inverse()", t1, t2, iter);
	//**********************************


	// Different conversions
	Point ppp(3,2,1);
	Vector nnn;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			nnn = ppp;
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Normal -> Point", t1, t2, iter);
	//**********************************

	std:: cout << "T"<< std::endl;

	//**********************************
	iter = 0;
	t1 = get_timestamp();
	while (iter < 10 ) {
		end = 0;
		while (end < 1000000) {
			nnn = Vector(ppp);
			++end;
		}
		++iter;
	}
	t2 = get_timestamp();
	print_time("Normal(Point) constructor", t1, t2, iter);
	//**********************************

	std::cout << std::endl << std::endl;
	std:: cout << "All OK!!!"<< std::endl << std::endl;

	return 0;
}           
