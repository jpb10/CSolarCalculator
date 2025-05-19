/*
 * SolarCalculator Library for Arduino (C version)
 *
 * This library provides functions to calculate the times of sunrise, sunset,
 * solar noon, twilight (dawn and dusk), Sun's apparent position in the sky,
 * equation of time, etc.
 *
 * Most formulae are taken from Astronomical Algorithms by Jean Meeus and
 * optimized for Arduino.
 */

#ifndef CSOLARCALCULATOR_H
#define CSOLARCALCULATOR_H

#ifdef __cplusplus
extern "C" {
#endif

#define SUNRISESET_ALTITUDE (-0.8333)
#define CIVIL_DAWNDUSK_ALTITUDE (-6.0)
#define NAUTICAL_DAWNDUSK_ALTITUDE (-12.0)
#define ASTRONOMICAL_DAWNDUSK_ALTITUDE (-18.0)

struct JulianDay
{
    double JD; /* Julian day at 0h UT (JD ending in .5) */
    double m;  /* fractional day (decimal number between 0 and 1) */
};

typedef struct JulianDay JulianDay;

/*
 * Get Julian Day, given:
 *
 *  - Unix timestamp, elapsed seconds since 1970 Jan 01, 00:00:00 UTC, e.g.
 *      ts      946684800       2000 Jan 01, 00:00:00 UTC
 *
 *  - Calendar date and time (UTC), with valid ranges:
 *      year    [1901..2099]
 *      month   [1..12]
 *      day     [1..31]
 *      hour    [0..23]
 *      minute  [0..59]
 *      second  [0..60]
 */

JulianDay jdFromUnix(unsigned long unix_ts);
JulianDay jdFromDate(int year, int month, int day,
                     int hour, int minute, int second);

/*
 * Intermediate calculations
 *
 * NB. Time T is measured in Julian centuries
 *     (36525 ephemeris days from the epoch J2000.0)
 */

/* Utilities */
double wrapTo360(double angle);
double wrapTo180(double angle);

/* Julian centuries */
double calcJulianCent(JulianDay jd);

/* Solar coordinates */
double calcGeomMeanLongSun(double T);
double calcGeomMeanAnomalySun(double T);
double calcSunEqOfCenter(double T);
double calcSunRadVector(double T);
double calcMeanObliquityOfEcliptic(double T);
void calcSolarCoordinates(double T,
                          double *ra, double *dec);

/* Sidereal time at Greenwich */
double calcGrMeanSiderealTime(JulianDay jd);

/* Sun's position in the sky */
void equatorial2horizontal(double H, double dec, double lat,
                           double *az, double *el);
double calcHourAngleRiseSet(double dec, double lat, double h0);
double calcRefraction(double elev);

/*
 * Solar calculator
 */

/* Equation of time, in minutes of time */
void calcEquationOfTime(JulianDay jd,
                        double *E);

/* Sun's (geocentric) equatorial coordinates, in degrees and AUs */
void calcEquatorialCoordinates(JulianDay jd,
                               double *rt_ascension, double *declination,
                               double *radius_vector);

/* Sun's (topocentric) horizontal coordinates, in degrees */
void calcHorizontalCoordinates(JulianDay jd, double latitude, double longitude,
                               double *azimuth, double *elevation);

/* Find the times of sunrise, transit, and sunset, in hours */
void calcSunriseSunset(JulianDay jd, double latitude, double longitude,
                       double *transit, double *sunrise, double *sunset,
                       double altitude, int iterations);

#ifdef __cplusplus
}
#endif

#endif /* CSOLARCALCULATOR_H */
