#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>
#include "sgp4unit.h"
#include "sgp4ext.h"
#include "SpaceObject.h"

#define SECOND_IN_JDAYS 1.157401129603386e-05; // 1 second expressed in Julian Days. JDAY/s
#define JDAY_IN_SECONDS 86400.46863810098; // 1 Julian Day expressed in seconds. s/JDAY

/* HELPER FUNCTIONS THAT HANDLE COMMAND LINE ARGUMENTS. */
char* getCmdOption(char ** begin, char ** end, const std::string & option){
	/* Get the value of the command line option that has been entered after the specified flag.
	@param begin, end - pointers to the beginning and the end of the char array of the command line arguments, i.e. argv, argv+argc, respectively.
	@param option - the flag that is to be found.
	@return - value corresponding to the given flag or 0 if nothing was supplied after the flag.
	*/
    char ** itr = std::find(begin, end, option);
    if(itr != end && ++itr != end){
        return *itr;
    };
    return 0;
};

bool cmdOptionExists(char** begin, char** end, const std::string& option){
	/* Check if a given command line option has been supplied.
	@param begin, end - pointers to the beginning and the end of the char array of the command line arguments, i.e. argv, argv+argc, respectively.
	@param option - the flag that is to be found.
	@return - true or false depending on whether the option has or has not been supplied to the program, respectively.
	*/
    return std::find(begin, end, option) != end;
};

std::map<std::string, SpaceObject> createObjectsFromTLEfile(const char* fileName, double defaultRadius_RB, double defaultRadius_PL, double defaultRadius_DEB, double defaultRadius_Other, int wgs, const char * radiusFileName){
	/* Read a file that contains multiple three line elements and create SpaceObjects from them. Return them to a std::map where the keys will be objects' NORAD IDs (SSCs).
	@param fileName - name of the file from which to read the TLE data.
	@param defaultR_RB - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains R/B.
	@param defaultR_PL - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains P/L.
	@param defaultR_DEB - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element contains DEB.
	@param defaultR_Other - default object to be used for object that have no entry in the radiusFileName in metres and the 1st line of the three-line element does not contain R/B, P/L or DEB.
	@param wgs - wgs ellipsoid to use, 72, 84, or 721 for old WGS 72
	@param radiusFileName - name of the file from which to read the hard body radii of objects.
	@return - map with SpaceObjects and their SSCs as keys.
	*/
	std::map<std::string, SpaceObject> SpaceObjects = std::map<std::string, SpaceObject>(); // A map containing SpaceObjects with their SSCs as keys.

	std::ifstream TLEfileStream(fileName); // File streams that are necessary.
	std::ifstream radiusFileStream(radiusFileName);
	
	/* Read the hard body radii to assign to objects when they are being created. */
	std::map<std::string, double>* hardBodyRadiiPtr = new std::map<std::string, double>();
	std::string radiusLine; // Currently read line with object's NORAD ID and hard body radius in m.

	int counterRadii = 0; // Current line of the TLE, 1 or 2.
	while( std::getline(radiusFileStream, radiusLine) ){
		std::istringstream radStringStream(radiusLine);
		if(counterRadii>5){ // Skip the header.
			int NORAD_ID; double radius; std::string NORAD_ID_STR;
			if (!(radStringStream >> NORAD_ID >> radius)) { break; } // Read the NORAD ID and the corresponding radius.
			
			std::stringstream strstream; // Convert the NORAD ID to a string.
			strstream << NORAD_ID; strstream>>NORAD_ID_STR;
			
			hardBodyRadiiPtr->insert ( std::pair<std::string,double>( NORAD_ID_STR, radius ) );
		}else{
			counterRadii+=1;
		};
	}

	/* Read the TLEs and create the SpaceObjects. */
	std::string TLEline; // Currently read TLE line.
	std::vector<std::string> currentTLE; currentTLE.resize(3); // Assembled TLE (second and third lines) and the object name (first line).

	int counterTLEs = 0; // Current line of the TLE, 0, 1 or 2.
	while( std::getline(TLEfileStream, TLEline) ){
		std::istringstream iss(TLEline);
		currentTLE.at(counterTLEs) = TLEline; // Add a new line.
		counterTLEs+=1;
		if(counterTLEs == 3){ // Read a whole TLE.
			counterTLEs = 0;
			SpaceObject tempSO = SpaceObject( currentTLE, 1, 1.0, wgs ); // Here one can change how covariance for the object is set (important for conjunctions) and whether to simulate changed atmospheric density (lieanrly proportional to B*).
			SpaceObjects.insert( std::pair<std::string,SpaceObject>( tempSO.NORAD_ID, tempSO) );
			try{
				SpaceObjects.at(tempSO.NORAD_ID).SetHardBodyRadius( hardBodyRadiiPtr->at(tempSO.NORAD_ID) ); // Assign the radius value.
				//std::cout<<"Created an object for "<<tempSO.NORAD_ID<<" with radius of "<<objectsPtr->at(tempSO.NORAD_ID).hardBodyRadius<<" km."<<std::endl;
			}catch(...){ // Radius for this object does not exist in the database, use a default value for this object type.
				if(currentTLE.at(0).find("R/B") != std::string::npos){
					SpaceObjects.at(tempSO.NORAD_ID).SetHardBodyRadius( defaultRadius_RB );
				}else if(currentTLE.at(0).find("P/L") != std::string::npos){
					SpaceObjects.at(tempSO.NORAD_ID).SetHardBodyRadius( defaultRadius_PL );
				}else if(currentTLE.at(0).find("DEB") != std::string::npos){
					SpaceObjects.at(tempSO.NORAD_ID).SetHardBodyRadius( defaultRadius_DEB );
				}else{
					SpaceObjects.at(tempSO.NORAD_ID).SetHardBodyRadius( defaultRadius_Other );
				};
			};
		};
	}

	char buff[25]; time_t now = time (0);
    #ifdef _MSC_VER 
		struct tm timeinfo;
		localtime_s(&timeinfo, &now);
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
	#else
		strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
	#endif
	std::cout<<buff<<": Generated "<<SpaceObjects.size()<<" objects."<<std::endl;
	
	return SpaceObjects;
};

/* MAIN PROGRAM THAT HANDLES INITALISATION OF EVERYTHING AS WELL AS MANAGES CONJUNCTION DETECTION AND RECORDING. */
int main(int argc, char *argv[]){

	std::string VERSION = "1.2.1";
	/* Check whether to display help. Terminate the program afterwards. */
	if( cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help") ){
		std::cout<<"Help for Conjunction Detection and Assessment v. "<<VERSION<<std::endl;
		std::cout<<"Argument flag / expected value type following the flag / default value / meaning"<<std::endl
			<<"-h || --help / NONE / NONE / Display this help information and terminate the program."<<std::endl
			<<"-TLE / const char * / spaceTrack3LEs_23_10_2013.txt / File from which the three line elements will be read and used in the analysis."<<std::endl
			<<"-radius / const char * / stkRadiusFile.rad / File from which the objects' radii will be read. Shall contain two columns of SSCs and corresponding radii in metres."<<std::endl
			<<"-wgs / integer 72, 84, or 721 for old WGS 72 / 84 / Change the WGS Earth ellipsoid type used."<<std::endl
			<<"-points / integer / 100 / Number of equispaced epochs at which to propagate the objects between jdayStart and jdayStop."<<std::endl
			<<"-rRB / double / 1.7691 / Default radius for the objects marked as R/B (rocket bodies) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object."<<std::endl
			<<"-rPL / double / 1.0350 / Default radius for the objects marked as P/L (payloads) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object."<<std::endl
			<<"-rDEB / double / 0.1558 / Default radius for the objects marked as DEB (debris objects) in the three-line element file. Default value comes from MASTER 2009 population average for this type of object."<<std::endl
			<<"-rOther / double / 0.3470 / Default radius for the objects not marked as R/B, P/L or DEB in the three-line element file. Default value comes from MASTER 2009 population average for this type of object."<<std::endl
			<<"-jdayStart / double / 23-10-2013 3:58:20.413924 / Julian day from which to start the analysis."<<std::endl
			<<"-jdayStop / double / 24-10-2013 3:58:20.413924 / Julian day at which to stop the analysis."<<std::endl;
		return 0;
	};

	/* Go through the command line arguments and set the supplied values or the default ones. */
	const char* TLEfile; // Objects' file.
	if( cmdOptionExists(argv, argv+argc, "-TLE") ){
		TLEfile = getCmdOption(argv, argv + argc, "-TLE");
		std::cout<<"Using object TLEs from "<<TLEfile<<std::endl;
	}else{
		TLEfile = "spaceTrack3LEs_23_10_2013.txt";
	};

	const char* RadFile; // Objects' radii file.
	if( cmdOptionExists(argv, argv+argc, "-radius") ){
		RadFile = getCmdOption(argv, argv + argc, "-radius");
		std::cout<<"Using object radii from "<<RadFile<<std::endl;
	}else{
		RadFile = "stkRadiusFile.rad";
	};

	int wgs; // WGS ellipsoid.
	if( cmdOptionExists(argv, argv+argc, "-wgs") ){
		std::basic_istringstream<char> wgsSS( getCmdOption(argv, argv + argc, "-wgs") );
		wgsSS >> wgs;
	}else{
		wgs=84;
	};

	int NoPoints; // Number of points at which to propagate the objects.
	if( cmdOptionExists(argv, argv+argc, "-points") ){
		std::basic_istringstream<char> pointsSS( getCmdOption(argv, argv + argc, "-points") );
		pointsSS >> NoPoints;
	}else{
		NoPoints=100;
	};

		double rRB, rPL, rDEB, rOther; // Default object radii.
	if( cmdOptionExists(argv, argv+argc, "-rRB") ){
		std::basic_istringstream<char> radRB_SS( getCmdOption(argv, argv + argc, "-rRB") );
		radRB_SS >> rRB;
	}else{
		rRB = 1.7691;
	};
	if( cmdOptionExists(argv, argv+argc, "-rPL") ){
		std::basic_istringstream<char> radPL_SS( getCmdOption(argv, argv + argc, "-rPL") );
		radPL_SS >> rPL;
	}else{
		rPL = 1.0350;
	};
	if( cmdOptionExists(argv, argv+argc, "-rDEB") ){
		std::basic_istringstream<char> radDEB_SS( getCmdOption(argv, argv + argc, "-rDEB") );
		radDEB_SS >> rDEB;
	}else{
		rDEB = 0.1558;
	};
	if( cmdOptionExists(argv, argv+argc, "-rOther") ){
		std::basic_istringstream<char> radOther_SS( getCmdOption(argv, argv + argc, "-rOther") );
		radOther_SS >> rOther;
	}else{
		rOther = 0.3470;
	};

		double ANALYSIS_INTERVAL_START; // Analysis start
	if( cmdOptionExists(argv, argv+argc, "-jdayStart") ){
		std::basic_istringstream<char> startSS( getCmdOption(argv, argv + argc, "-jdayStart") );
		startSS >> ANALYSIS_INTERVAL_START;
	}else{
		jday(2013,10,23,3,58,20.413924, ANALYSIS_INTERVAL_START);
	};

	double ANALYSIS_INTERVAL_STOP; // Analysis stop
	if( cmdOptionExists(argv, argv+argc, "-jdayStop") ){
		std::basic_istringstream<char> stopSS( getCmdOption(argv, argv + argc, "-jdayStop") );
		stopSS >> ANALYSIS_INTERVAL_STOP;
	}else{
		jday(2013,10,24,3,58,20.413924, ANALYSIS_INTERVAL_STOP);
	};

	/* Run the actual program. */
	std::map<std::string, SpaceObject> SpaceObjects = createObjectsFromTLEfile(TLEfile, rRB, rPL, rDEB, rOther, wgs, RadFile); // Read the objects' TLEs from the file.

	std::vector<double> AnalysisJDAYs = std::vector<double>(NoPoints, 0.0); // Epochs at which the objects will be propagated.
	linspace(&AnalysisJDAYs, ANALYSIS_INTERVAL_START, ANALYSIS_INTERVAL_STOP, &NoPoints);

	std::vector<double> position = std::vector<double>(3, 0.0); // Cartesian position of a given object in TEME of Date in km.
	std::vector<double> velocity = std::vector<double>(3, 0.0); // Cartesian velocity of a given object in TEME of Date in km/sec.

	for(std::map<std::string, SpaceObject>::iterator outerIt = SpaceObjects.begin(); outerIt != SpaceObjects.end(); ++outerIt){
		for(std::vector<int>::size_type i=0; i<AnalysisJDAYs.size(); i++){
			outerIt->second.PropagateJDAY(&position, &velocity, AnalysisJDAYs.at(i));
		};
	};
	
	return 0;
}
