#include "SpaceObject.h"

SpaceObject::SpaceObject(void){
	/* Default constructor, does nothing special. */
};

SpaceObject::SpaceObject(std::vector<std::string> threeLineElement, int COV, double bstarMultiplier, int wgsEllipsoid, double objectRadius, char opsMode){
	/* Create an object that represents something orbiting the Earth - a convenience interface for SGP4 propagator by D. Vallado.
	@param threeLineElement - a vector of strings that cointains the object name (line 0) and the first and second lines of the TLE file for a given object. Can have trailing end of line characters and must have the leading 1 and 2.
	@param COV - type of covariance to be used for the objects. If 1 an uncertainty sphere 1:1:1 will be used and only maximum collision probability will be computed.
		If 2 a covariance will be estimated.
	@param bstarMultiplier - a factor by which the B* coefficient of a satellite will be multiplied by, used to simulate atmospheric density changes.
	@param wgsEllispoid - specifies which WGS Earth ellipsoid to use, possible options are 72, 84 and 721 (old WGS72).
	@param objectRadius - radius of the object in metres, used to compute the probability of collisions.
	@param opsMode - operation mode of the SGP4, either i (improved) or a (afspc).
	*/
	char OPSMODE = opsMode; // Improved operation mode of SGP4, should result in a smoother behaviour.
	earthEllipsoid = wgsEllipsoid;
	BstarMultiplier = bstarMultiplier;

	/* Initialise SGP4-specific attributes. */
	revnum=0; elnum=0;
	
	// These are defined in SGP4EXT, will need to uncomment them when not using it.
	//const double deg2rad = PI/180.0;         // Conversion factor. 0.0174532925199433
	//const double xpdotp = 1440.0/(2.0*PI);  // Minutes in a dayper a full revolution = 229.1831180523293. Unit: min/rad
	sgp4Sat = elsetrec(); // SGP4 satellite struct that is used for propagation using the SGP4.
	
	/* Define the desired EGS ellipsoid. */
	if(wgsEllipsoid==84){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==72){
		gravityModel = wgs84;
	}else if(wgsEllipsoid==721){
		gravityModel = wgs72;
	};
	getgravconst( gravityModel, tumin, mu, Re, xke, j2, j3, j4, j3oj2 ); // Get values for these variables.

	/* Convert the TLE into the expected format. */
	char TLE1[130];	char TLE2[130]; // TLE lines - char arrays of the appropriate length that are expected by sgp4init.

	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		strcpy_s(TLE1, threeLineElement.at(1).c_str());
		strcpy_s(TLE2, threeLineElement.at(2).c_str());
	#else
		std::strcpy(TLE1, threeLineElement.at(1).c_str());
		std::strcpy(TLE2, threeLineElement.at(2).c_str());
	#endif
	
	/* Set the implied decimal points since doing a formated read
     * fixes for bad input data values (missing, ...). */
    for (int j = 10; j <= 15; j++){
        if (TLE1[j] == ' ')
			TLE1[j] = '_';
	};
	if (TLE1[44] != ' ')
		TLE1[43] = TLE1[44];
	TLE1[44] = '.';
	if (TLE1[7] == ' ')
		TLE1[7] = 'U';
    if (TLE1[9] == ' ')
		TLE1[9] = '.';
	for (int j = 45; j <= 49; j++){
		if (TLE1[j] == ' ')
			TLE1[j] = '0';
	};
    if (TLE1[51] == ' ')
		TLE1[51] = '0';
    if (TLE1[53] != ' ')
		TLE1[52] = TLE1[53];
	TLE1[53] = '.';
    TLE2[25] = '.';
    for (int j = 26; j <= 32; j++){
        if (TLE2[j] == ' ')
			TLE2[j] = '0';
	};
    if (TLE1[62] == ' ')
        TLE1[62] = '0';
    if (TLE1[68] == ' ')
        TLE1[68] = '0';

	/* Parse the TLE. */
	#ifdef _MSC_VER // Depending on the compiler being used utilise the appropriate function to copy the TLE lines.
		sscanf_s(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification, sizeof(char),intldesg, 11*sizeof(char),&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf_s(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#else
		sscanf(TLE1,"%2d %5ld %1c %10s %2d %12lf %11lf %7lf %2d %7lf %2d %2d %6ld ",&cardnumb,&sgp4Sat.satnum,&classification,intldesg,&sgp4Sat.epochyr,&sgp4Sat.epochdays,&sgp4Sat.ndot,&sgp4Sat.nddot,&nexp,&sgp4Sat.bstar,&ibexp,&numb,&elnum);
		sscanf(TLE2,"%2d %5ld %9lf %9lf %8lf %9lf %9lf %11lf %6ld %lf %lf %ld \n",&cardnumb,&sgp4Sat.satnum,&sgp4Sat.inclo,&sgp4Sat.nodeo,&sgp4Sat.ecco,&sgp4Sat.argpo,&sgp4Sat.mo,&sgp4Sat.no,&revnum);
	#endif
	
	std::stringstream strstream; // Record the NORAD ID as a string as well.
	strstream << sgp4Sat.satnum; strstream >> NORAD_ID;
		
	/* Iniitalise the SGP4 propagator variables and the propagator itself for this object. */
	// Find the TLE epoch in Julian days.
    int year = 2000 + sgp4Sat.epochyr; // N.B. this won't work for historic TLEs from before 2000.

	int epMon, epDay, epHr, epMin; double epSec; // TLE epoch components.
	days2mdhms(year, sgp4Sat.epochdays, epMon, epDay, epHr, epMin, epSec);
	jday(year, epMon, epDay, epHr, epMin, epSec, sgp4Sat.jdsatepoch); 

	// Find no, ndot, nddot.
    sgp4Sat.no   = sgp4Sat.no / xpdotp; // When using SGP4EXT this will already be in correct units, anyway multiply by: rad/min
    sgp4Sat.nddot= sgp4Sat.nddot * pow(10.0, nexp);
    sgp4Sat.bstar= sgp4Sat.bstar * pow(10.0, ibexp) * bstarMultiplier; // Multiply by the factor that allows variations in solar activity to be synthesised.
    
	// Convert to sgp4 units.
    sgp4Sat.a    = pow( sgp4Sat.no*tumin , (-2.0/3.0) );
    sgp4Sat.ndot = sgp4Sat.ndot  / (xpdotp*1440.0);  //* ? * minperday
    sgp4Sat.nddot= sgp4Sat.nddot / (xpdotp*1440.0*1440);
    
	// Find standard orbital elements.
    sgp4Sat.inclo = sgp4Sat.inclo  * deg2rad;
    sgp4Sat.nodeo = sgp4Sat.nodeo  * deg2rad;
    sgp4Sat.argpo = sgp4Sat.argpo  * deg2rad;
    sgp4Sat.mo    = sgp4Sat.mo     * deg2rad;

	// Iniitlaise the SGP4 for this TLE.
	sgp4init( gravityModel, OPSMODE, sgp4Sat.satnum, sgp4Sat.jdsatepoch-2433281.5, sgp4Sat.bstar, sgp4Sat.ecco, sgp4Sat.argpo, sgp4Sat.inclo, sgp4Sat.mo, sgp4Sat.no, sgp4Sat.nodeo, sgp4Sat);

	/* Call the propagator to get the initial state. */
	double r0[3]; double v0[3]; // Create arrays that SGP4 expects
    sgp4(gravityModel, sgp4Sat,  0.0, r0,  v0); // Third argument is time since epoch in minutes - 0.0 gives the initial position and velocity.
	currentPos.resize( sizeof(r0)/sizeof(r0[0]) ); currentVelo.resize( sizeof(v0)/sizeof(v0[0]) ); // Resize the vectors to be of appropriate length.
	currentPos = std::vector<double>(r0, r0+3); currentVelo = std::vector<double>(v0, v0+3); // Record the values of the initial position and velocity in the vectors.
	currentEpochJDAY = sgp4Sat.jdsatepoch;
	TLEepochJDAY = sgp4Sat.jdsatepoch; // In Julian Days.

	/* Get other values of interest, mianly orbital elements. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(r0, v0, mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	SemiMajorAxis = SMA; Eccentricity = ECC; Inclination = INCL; MeanAnomaly = MA; LongAscendingNode = LAN; ArgumentOfPerigee = ARGP;
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.

	hardBodyRadius = objectRadius/1000.; // Default value in km.

	if(COV==2){
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));
	}else if(COV==1){
		FullCovarianceMatrixRTC = std::vector< std::vector<double> >(6, std::vector<double>(6, 0.0)); // Initialise the covariance matrices.
		PositionCovarianceMatrixRTC = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0));

		// Initialise the diagonal entries of the covariance matrices to simulate the desired aspect ratio.
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0; FullCovarianceMatrixRTC.at(3).at(3)=1.0; FullCovarianceMatrixRTC.at(4).at(4)=1.0; FullCovarianceMatrixRTC.at(5).at(5)=1.0;
		FullCovarianceMatrixRTC.at(0).at(0)=1.0; FullCovarianceMatrixRTC.at(1).at(1)=1.0; FullCovarianceMatrixRTC.at(2).at(2)=1.0;
	}
};

void SpaceObject::PropagateJDAY(std::vector<double>* posPtr, std::vector<double>* veloPtr, double JDAY, bool updateCurrentState){
	/* Propagate the satellite to the specified Julain day and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param JDAY - Julian day at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	*/
	double r[3]; double v[3]; // Create arrays that SGP4 expects
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0;
    sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);
	
	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length here as well - they are Cartesian components.
			currentPos.at(i) = r[i];
			currentVelo.at(i) = v[i];
		};
	};
};

void SpaceObject::Propagate(std::vector<double>* posPtr, std::vector<double>* veloPtr, int year, int month, int day, int hour, int minute, double second, bool updateCurrentState){
	/* Propagate the satellite to the specified Universal Time epocj and output its position and velocity at that epoch to currentPosPtr
	and currentVeloPtr.
	@param posPtr, veloPtr - pointers to std::vector<double> of length 3 that contain current positions and velocities of the object.
	@param year, month, day, hour, minute, second - Universal Time at which the object's state is to be computed.
	@param updateCurrentState - whether to update currentPos, currentVelo and currentEpochJDAY with the new values.
	*/	
	double JDAY = currentEpochJDAY; // Initialise with the last new location.
	jday(year, month, day, hour, minute, second, JDAY); // Convert desired epoch to Julian days.
	double minutesSinceEpoch = (JDAY-TLEepochJDAY)*1440.0; // And minutes since TLE epoch.
	
	double r[3]; double v[3]; // Propagate using arrays that SGP4 expects
	sgp4(gravityModel, sgp4Sat,  minutesSinceEpoch, r,  v);

	for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Record the position and velocity in the vectors.
			posPtr->at(i) = r[i];
			veloPtr->at(i) = v[i];
		};

	if(updateCurrentState){ // Update the state variables if desired.
		currentEpochJDAY = JDAY;
		for(std::vector<int>::size_type i = 0; i < currentPos.size(); i++){ // Both vectors have the same length - they are Cartesian components.
			currentPos[i] = r[i];
			currentVelo[i] = v[i];
		};
	};
};

void SpaceObject::CalculatePerigeeRadius(void){
	/* Update the currentPerigeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_PERIGEE_RADIUS = SMA*(1.0 - ECC); // In km.
};

void SpaceObject::CalculateApogeeRadius(void){
	/* Update the currentApogeeRadius based on currentPos and currentVelo values. */
	double SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper;
	rv2coe(&currentPos[0], &currentVelo[0], mu, SLR, SMA, ECC, INCL, LAN, ARGP, TA, MA, arglat, truelon, lonper );
	CURRENT_APOGEE_RADIUS = SMA*(1.0 + ECC); // In km.
};

void SpaceObject::SetHardBodyRadius(double bodyRadius){
	/* Set the hardBodyRadius attribute to the specified value in m. */
	hardBodyRadius = bodyRadius/1000.; // Convert to km.
};

std::vector< std::vector<double> > SpaceObject::InertialToRTC(std::vector<double>* positionPtr, std::vector<double>* velocityPtr){
    /* Compute a rotation matrix from an inertial frame of reference (TEME, ECI or similar) to 
    Radial - Transverse - In-track (also called Radial - Along-track - Cross-ctrack) frame.
    @param position, velcity - dimensional position and velocity in the inertial frame. Must have the same units, assumed km and km/sec.
    @return - 3x3 rotation matrix from inertial to RTC frames - DOES NOT INCLUDE TRANSLATION BETWEEN ORIGINS.
	*/
	/* Base vectors of the RTC frame expressed in the inertial frame. */
    std::vector<double> radialUnit_inertial = vectorMultiplyByScalar( positionPtr, 1.0/vectorMagnitude(positionPtr) ); // Unit radial vector expressed in inertial frame.
	
    std::vector<double> crossUnit_inertial = std::vector<double>(3, 0.0); // Cross-track unit vector expressed in the inertial reference frame.
	crossProduct( &crossUnit_inertial, positionPtr, velocityPtr);
	unitVector( &crossUnit_inertial ); // Make unit.
    
	std::vector<double> transverseUnit_inertial = std::vector<double>(3, 0.0); // Transverse unit vector expressed in the inertial reference frame.
	crossProduct( &transverseUnit_inertial, &radialUnit_inertial, &crossUnit_inertial);
	unitVector( &transverseUnit_inertial ); // Make unit.

	/* Base vectors of the inertial frame expressed in the inertial frame. */
    std::vector<double> xAxisUnitVectorInertial = std::vector<double>(3, 0.0); xAxisUnitVectorInertial.at(0)=1.0;
    std::vector<double> yAxisUnitVectorInertial = std::vector<double>(3, 0.0); yAxisUnitVectorInertial.at(1)=1.0;
    std::vector<double> zAxisUnitVectorInertial = std::vector<double>(3, 0.0); zAxisUnitVectorInertial.at(2)=1.0;

	/* Formulate the rotation matrix from the inertial to RTC frame. */
    std::vector< std::vector<double> > rotationMatrix = std::vector< std::vector<double> >(3, std::vector<double>(3, 0.0) ); // Full 3x3 rotation matrix from inertial to RTC coordinate frames.
	rotationMatrix.at(0).at(0) = dotProduct( &radialUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(0).at(1) = dotProduct( &radialUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(0).at(2) = dotProduct( &radialUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(1).at(0) = dotProduct( &transverseUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(1).at(1) = dotProduct( &transverseUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(1).at(2) = dotProduct( &transverseUnit_inertial, &zAxisUnitVectorInertial );

	rotationMatrix.at(2).at(0) = dotProduct( &crossUnit_inertial, &xAxisUnitVectorInertial );
	rotationMatrix.at(2).at(1) = dotProduct( &crossUnit_inertial, &yAxisUnitVectorInertial );
	rotationMatrix.at(2).at(2) = dotProduct( &crossUnit_inertial, &zAxisUnitVectorInertial );

	return rotationMatrix;
};

SpaceObject::~SpaceObject(void){
 /* Deconstructor, does nothing special.*/
};
