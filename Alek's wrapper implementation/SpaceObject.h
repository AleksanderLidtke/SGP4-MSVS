/*
A class that provides a convenient interface to the SGP4 propagator.
Based on code by David Vallado freely available on clestrak.com (last accessed 03/02/2014).

@version 1.1.1
@since 24/07/2014 10:29:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:
23/07/2014 - 1.1.0 - Added the covariance parameter that enables a fixed covariance matrix to be set.
24/07/2014 - 1.1.1 - fixed a bug where the NORAD_ID attribute would be changed hence causing the entire object creation framework to fail.
*/
#pragma once

#include <iostream>
#include <stdio.h>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "sgp4ext.h"
#include "sgp4unit.h"
#include "VectorOperations.h"

#define PI 3.14159265358979323846

class SpaceObject{
private:
	std::vector< std::vector<double> > InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr);
public:
	/* Methods. */
	SpaceObject(void);
	SpaceObject(std::vector<std::string> threeLineElement, int COV=1, double bstarMultiplier=1.0, int wgsEllipsoid=84, double objectRadius=5.0, char opsMode='i');
	~SpaceObject(void);

	void PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* vloPtr, double JDAY, bool updateCurrentState=false);
	void Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second,  bool updateCurrentState=false);
	void CalculatePerigeeRadius(void);
	void CalculateApogeeRadius(void);
	void SetHardBodyRadius(double bodyRadius);

	/* Propagation-sepcific attrributes. */
	char OPSMODE; // Improved operation mode of SGP4, should result in a smoother behaviour.
	double tumin, mu, Re, xke, j2, j3, j4, j3oj2; // Gravity model properties.
	gravconsttype gravityModel; // Gravity model being used - WGS84, WGS72 or old WGS72.
	elsetrec sgp4Sat; // SGP4 propagation object.
	char classification; // Object classification (all here will be U - unclassified.
	char intldesg[11]; // International designator - launch year, launch number, launch piece.
	double BstarMultiplier; // Multiplication factor of the TLEs' B*, used to simulate variations in solar activity.

	/* SGP4-specific attributes. */
	int cardnumb, numb;
    long revnum, elnum;
	int nexp, ibexp;

	/* General orbital elements and other attributes - may be updated as propagation progresses but initialised with the values at TLE epoch. */
	std::vector<double> currentPos, currentVelo; // Current position and velocity.
	double CURRENT_PERIGEE_RADIUS, CURRENT_APOGEE_RADIUS, currentEpochJDAY; // Corresponding radii and epoch.
	double hardBodyRadius, TLEepochJDAY; // Physical radius of the object (used to compute the collision probability) and epoch of the TLE that was used to create it.
	int earthEllipsoid; // Year of the WGS ellipsoid used to propgate the object, 72 or 84.
	std::string NORAD_ID; // SSC/NORAD ID number of the object.
	double SemiMajorAxis, Eccentricity, Inclination, MeanAnomaly, LongAscendingNode, ArgumentOfPerigee; // Classical orbital elements at the epoch of the TLE.

	/* Covariance matrix related attributes. */
	std::vector< std::vector<double> > FullCovarianceMatrixRTC; // 6x6 position and velocity covariance matrix in km and km/sec.
	std::vector< std::vector<double> > PositionCovarianceMatrixRTC; // 3x3 position covariance matrix in km.
};
