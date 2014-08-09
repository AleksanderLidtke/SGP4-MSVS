#ifndef _EopSpw_h_
#define _EopSpw_h_
/* ---------------------------------------------------------------------
*
*                              EopSpw.h
*
*  this file contains routines to read the eop and space weather data
*  from the cssi files on celestrak. it also includes routines to
*  interpolate and convert between the values.
*
*  current :
*             8 aug 14 alek lidtke
*                           removed hard-coded include paths
*  changes :
*             7 aug 14  david vallado
*                           convert to msvs2013 c++
*            14 dec 05  david vallado
*                           misc fixes
*            21 oct 05  david vallado
*                           original version
*       ----------------------------------------------------------------      */

#pragma once

#include <math.h>
#include <io.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "astMath.h"
#include "ast2body.h"
#include "astTime.h"
//#include "msiscom.h"
//#include "msis00.h"


//using namespace std;
//
//using namespace System;
//

/*    *****************************************************************
*     type definitions
*     *****************************************************************     */

// note the deltapsi/deltaeps are for iau76/fk5 from eop
// x/y/s/a and deltapsi/deltaeps are for computational speed in nutation calcs

typedef struct eopdata
  {
       double xp,   yp,  dut1, lod,  ddpsi,    ddeps,    dx,   dy;
       int    year, mon, day,  mjd,  dat;
       double x,    y,   s,    deltapsi, deltaeps;
  } eopdata;

typedef struct spwdata
  {
       double adjf10,  adjctrf81, adjlstf81, obsf10,   obsctrf81, obslstf81, cp;
       int    year,    mon,       day,       bsrn,     nd,        avgap,     c9,
              isn,     q,         aparr[8],  kparr[8], sumkp;
  } spwdata;

const int eopsize = 2500; // 25000 if from 62
const int spwsize = 2500; // 25000 if from 62


namespace EopSpw {

	//	public ref class EopSpwCl
	//	{

	// make sure they are all visible
	//	public:


	void initspw
		(
		spwdata spwarr[spwsize],
		double& jdspwstart
		);

	void initeop
		(
		eopdata eoparr[eopsize],
		double& jdeopstart
		);

	void findeopparam
		(
		double  jd, double mfme, char interp,
		eopdata eoparr[eopsize], double jdeopstart,
		double& dut1, int& dat,
		double& lod, double& xp, double& yp,
		double& ddpsi, double& ddeps, double& dx, double& dy,
		double& x, double& y, double& s,
		double& deltapsi, double& deltaeps
		);

	void findatmosparam
		(
		double jd, double mfme, char interp, char fluxtype, char f81type, char inputtype,
		spwdata spwarr[spwsize], double jdspwstart,
		double& f107, double& f107bar,
		double& ap, double& avgap, double aparr[8],
		double& kp, double& sumkp, double kparr[8]
		);

	double kp2ap
		(
		double kpin
		);

	double ap2kp
		(
		double apin
		);




	// this routine will probably move elsewhere, but fits here ok for now
	//void interfaceatmos
	//    (
	//       double jde, double mfme, double recef[3],
	//       char interp, char fluxtype, char f81type, char inputtype,
	//       msistype& msis00r,
	//       spwdata spwarr[spwsize], double jdspwstart
	//     );

	//	}; // class

};  // namespace

#endif
