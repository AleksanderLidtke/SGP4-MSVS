/*
A custom error class that will include a time stamp for the error.

@version 1.0.0
@since 20/04/2014 23:15:00
@author Aleksander Lidtke
@email al11g09@soton.ac.uk, alek_l@onet.eu

CHANGELOG:

*/
#pragma once

#include <time.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>

class FatalError : public std::runtime_error{
	public:
		FatalError(void):runtime_error("Fatal error"){}
		FatalError(std::string msg, std::string file, int line):runtime_error(msg.c_str()) {
			char buff[25]; time_t now = time (0);
			#ifdef _MSC_VER 
				struct tm timeinfo;
				localtime_s(&timeinfo, &now);
				strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", &timeinfo);
			#else
				strftime(buff, 25, "%Y-%m-%d %H:%M:%S.000", localtime (&now));
			#endif

			std::cout <<buff<<" ---> FATAL ERROR IN: " <<file<< " AT LINE " << line << ":\n\n"; }
};
