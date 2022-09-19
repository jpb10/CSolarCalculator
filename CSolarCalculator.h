//======================================================================================================================
// SolarCalculator Library for Arduino (C version)
//
// This library provides functions to calculate the times of sunrise, sunset, solar noon, twilight (dawn and dusk),
// Sun's apparent position in the sky, equation of time, etc.
//
// Most formulae are taken from Astronomical Algorithms by Jean Meeus and optimized for 8-bit AVR platform.
//======================================================================================================================

#ifndef CSOLARCALCULATOR_H
#define CSOLARCALCULATOR_H

#ifdef __cplusplus
extern "C" {
#endif

extern const double SUNRISESET_STD_ALTITUDE;
extern const double CIVIL_DAWNDUSK_STD_ALTITUDE;
extern const double NAUTICAL_DAWNDUSK_STD_ALTITUDE;
extern const double ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE;

struct JulianDay
{
    double JD;  // Julian day at 0h UT (JD ending in .5)
    double m;   // Fractional day, 0h to 24h (decimal number between 0 and 1)
};

struct JulianDay jdFromUnix(unsigned long utc);  // Unix time, i.e. seconds since 0h UT 1 January 1970
struct JulianDay jdFromDate(int year, int month, int day, int hour, int minute, int second);  // Calendar date (UTC)

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian centuries (36525 ephemeris days from the epoch J2000.0)
//======================================================================================================================

// Utilities
double wrapTo360(double angle);
double wrapTo180(double angle);

// Julian day
double fractionalDay(int hour, int minute, int second);
double calcJulianDay(int year, int month, int day);
double calcJulianCent(struct JulianDay jd);

// Solar coordinates
double calcGeomMeanLongSun(double T);
double calcGeomMeanAnomalySun(double T);
double calcSunEqOfCenter(double T);
double calcSunRadVector(double T);
double calcMeanObliquityOfEcliptic(double T);
void calcSolarCoordinates(double T, double *ra, double *dec);

// Sidereal time, Sun's position in the sky
double calcGrMeanSiderealTime(struct JulianDay jd);
void equatorial2horizontal(double H, double dec, double lat, double *az, double *el);
double calcRefraction(double elev);
double equationOfTimeSmart(double T);

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

// Equation of time, in minutes of time
void calcEquationOfTime(struct JulianDay jd, double *E);

// Sun's geocentric (as seen from the center of the Earth) equatorial coordinates, in degrees and AUs
void calcEquatorialCoordinates(struct JulianDay jd, double *rt_ascension, double *declination, double *radius_vector);

// Sun's topocentric (as seen from the observer's place on the Earth's surface) horizontal coordinates, in degrees
void calcHorizontalCoordinates(struct JulianDay jd, double latitude, double longitude,
                               double *azimuth, double *elevation);

// Find the times of sunrise, transit, and sunset, in hours
void calcSunriseSunset(struct JulianDay jd, double latitude, double longitude,
                       double *transit, double *sunrise, double *sunset,
                       double altitude, int iterations);

//======================================================================================================================
// Wrapper functions
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

void calcEquationOfTimeFromUnix(unsigned long utc, double *E);
void calcEquationOfTimeFromDate(int year, int month, int day, int hour, int minute, int second, double *E);

void calcEquatorialCoordinatesFromUnix(unsigned long utc,
                                       double *rt_ascension, double *declination, double *radius_vector);
void calcEquatorialCoordinatesFromDate(int year, int month, int day, int hour, int minute, int second,
                                       double *rt_ascension, double *declination, double *radius_vector);

void calcHorizontalCoordinatesFromUnix(unsigned long utc, double latitude, double longitude,
                                       double *azimuth, double *elevation);
void calcHorizontalCoordinatesFromDate(int year, int month, int day, int hour, int minute, int second,
                                       double latitude, double longitude, double *azimuth, double *elevation);

void calcSunriseSunsetFromUnix(unsigned long utc, double latitude, double longitude,
                               double *transit, double *sunrise, double *sunset);
void calcSunriseSunsetFromDate(int year, int month, int day, double latitude, double longitude,
                               double *transit, double *sunrise, double *sunset);

void calcCivilDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                               double *transit, double *dawn, double *dusk);
void calcCivilDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                               double *transit, double *dawn, double *dusk);

void calcNauticalDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                                  double *transit, double *dawn, double *dusk);
void calcNauticalDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                                  double *transit, double *dawn, double *dusk);

void calcAstronomicalDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                                      double *transit, double *dawn, double *dusk);
void calcAstronomicalDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                                      double *transit, double *dawn, double *dusk);

#ifdef __cplusplus
}
#endif

#endif  //CSOLARCALCULATOR_H
