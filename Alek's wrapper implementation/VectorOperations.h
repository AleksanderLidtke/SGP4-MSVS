/*
A set of functions to perform basic vector and matrix operations on 1D vectors or 2D matrices.

@version 1.0.0
@since 20/04/2014 23:15:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:

*/
#pragma once

#include <cmath>
#include <vector>
#include "FatalError.h"

void linspace(std::vector<double>* outputVecPtr, double start, double stop, int* noPoints);
void centralLinspace(std::vector<double>* outputVecPtr, double centre, double spacing, int* noPoints);
void vectorDifference(std::vector<double>* outputVecPtr, std::vector<double>* vect1Ptr, std::vector<double>* vec2Ptr);
double vectorMagnitudeSquared(std::vector<double>* vecPtr);
double vectorMagnitude(std::vector<double>* vecPtr);
void unitVector(std::vector<double>* vecPtr);
double dotProduct(std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr);
void crossProduct(std::vector<double>* outputVecPtr, std::vector<double>* vec1Ptr, std::vector<double>* vec2Ptr);
std::vector<double> vectorMultiplyByScalar(std::vector<double>* vecPtr, double coefficient);
std::vector< std::vector<double> > matrixMultiplyByScalar(std::vector< std::vector<double> >* matPtr, double coefficient);
std::vector <double> vectorMultiplyByMatrix(std::vector<double>* vecPtr, std::vector< std::vector<double> >* matPtr);
std::vector< std::vector<double> > matrixMultiplyByMatrix(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr);
std::vector< std::vector<double> > addMatrices(std::vector< std::vector<double> >* mat1Ptr, std::vector< std::vector<double> >* mat2Ptr);
double twoByTwoDeterminant( std::vector< std::vector<double> >* matPtr );
double threeByThreeDeterminant( std::vector< std::vector<double> >* matPtr );
std::vector< std::vector<double> > inverseTwoByTwo( std::vector< std::vector<double> >* matPtr );
std::vector< std::vector<double> > inverseThreeByThree( std::vector< std::vector<double> >* matPtr );
std::vector< std::vector<double> > transposeMatrix( std::vector< std::vector<double> >* originalPtr );