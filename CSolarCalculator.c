//======================================================================================================================
// SolarCalculator Library for Arduino (C version)
//
// This library provides functions to calculate the times of sunrise, sunset, solar noon, twilight (dawn and dusk),
// Sun's apparent position in the sky, equation of time, etc.
//
// Most formulae are taken from Astronomical Algorithms by Jean Meeus and optimized for 8-bit AVR platform.
//======================================================================================================================

#ifdef ARDUINO
#include <Arduino.h>
#else
#include <math.h>
#endif

#include "CSolarCalculator.h"

const double SUNRISESET_STD_ALTITUDE = -0.8333;
const double CIVIL_DAWNDUSK_STD_ALTITUDE = -6.0;
const double NAUTICAL_DAWNDUSK_STD_ALTITUDE = -12.0;
const double ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE = -18.0;

struct JulianDay jdFromUnix(unsigned long utc)
{
    struct JulianDay jd;
    jd.JD = (unsigned long)(utc / 86400) + 2440587.5;
    jd.m = (utc % 86400) / 86400.0;
    return jd;
}

struct JulianDay jdFromDate(int year, int month, int day, int hour, int minute, int second)
{
    struct JulianDay jd;
    jd.JD = calcJulianDay(year, month, day);
    jd.m = fractionalDay(hour, minute, second);
    return jd;
}

//======================================================================================================================
// Intermediate calculations
//
// Time T is measured in Julian centuries (36525 ephemeris days from the epoch J2000.0)
//======================================================================================================================

#ifndef ARDUINO
double radians(double deg)
{
    return deg * M_PI / 180;
}

double degrees(double rad)
{
    return rad * 180 / M_PI;
}
#endif

double wrapTo360(double angle)
{
    angle = fmod(angle, 360);
    if (angle < 0) angle += 360;
    return angle;  // [0, 360)
}

double wrapTo180(double angle)
{
    angle = wrapTo360(angle + 180);
    return angle - 180;  // [-180, 180)
}

double fractionalDay(int hour, int minute, int second)
{
    return (hour + minute / 60.0 + second / 3600.0) / 24;
}

// Valid from 1901 to 2099, Van Flandern & Pulkkinen (1979)
double calcJulianDay(int year, int month, int day)
{
    return 367.0 * year - (int)(7 * (year + (month + 9) / 12) / 4) + (int)(275 * month / 9) + day + 1721013.5;
}

double calcJulianCent(struct JulianDay jd)
{
    return ((jd.JD - 2451545) + jd.m) / 36525;
}

double calcGeomMeanLongSun(double T)
{
    return wrapTo360(280.46646 + T * 36000.76983);  // in degrees
}

double calcGeomMeanAnomalySun(double T)
{
    return wrapTo360(357.52911 + T * 35999.05029);  // in degrees
}

double calcSunEqOfCenter(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    return sin(radians(M)) * 1.914602 + sin(2 * radians(M)) * 0.019993;  // in degrees
}

double calcSunRadVector(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    return 1.00014 - 0.01671 * cos(radians(M)) - 0.00014 * cos(2 * radians(M));  // in AUs
}

double calcMeanObliquityOfEcliptic(double T)
{
    return 23.4392911 - T * 0.0130042;  // in degrees
}

// Mean geocentric equatorial coordinates, accurate to ~1 arcminute
void calcSolarCoordinates(double T, double *ra, double *dec)
{
    double L0 = calcGeomMeanLongSun(T);
    double C = calcSunEqOfCenter(T);
    double L = L0 + C - 0.00569;  // corrected for aberration

    double eps = calcMeanObliquityOfEcliptic(T);
    *ra = degrees(atan2(cos(radians(eps)) * sin(radians(L)), cos(radians(L))));  // [-180, 180)
    *dec = degrees(asin(sin(radians(eps)) * sin(radians(L))));
}

double calcGrMeanSiderealTime(struct JulianDay jd)
{
    double GMST = wrapTo360(100.46061837 + 0.98564736629 * (jd.JD - 2451545));
    return wrapTo360(GMST + 360.985647 * jd.m);  // in degrees
}

void equatorial2horizontal(double H, double dec, double lat, double *az, double *el)
{
    *az = degrees(atan2(sin(radians(H)), cos(radians(H)) * sin(radians(lat)) - tan(radians(dec)) * cos(radians(lat))));
    *el = degrees(asin(sin(radians(lat)) * sin(radians(dec)) + cos(radians(lat)) * cos(radians(dec)) * cos(radians(H))));
}

// Approximate atmospheric refraction correction, in degrees
double calcRefraction(double elev)
{
    if (elev < -0.575)
        return -20.774 / tan(radians(elev)) / 3600;  // Zimmerman (1981)
    else
        return 1.02 / tan(radians(elev + 10.3 / (elev + 5.11))) / 60;  // S??mundsson (1986)
}

//======================================================================================================================
// Solar calculator
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

// Equation of time, in minutes of time
void calcEquationOfTime(struct JulianDay jd, double *E)
{
    double T = calcJulianCent(jd);
    double L0 = calcGeomMeanLongSun(T);

    double ra, dec;
    calcSolarCoordinates(T, &ra, &dec);

    *E = 4 * wrapTo180(L0 - 0.00569 - ra);
}

// Sun's geocentric (as seen from the center of the Earth) equatorial coordinates, in degrees and AUs
void calcEquatorialCoordinates(struct JulianDay jd, double *rt_ascension, double *declination, double *radius_vector)
{
    double T = calcJulianCent(jd);
    calcSolarCoordinates(T, rt_ascension, declination);

    *rt_ascension = wrapTo360(*rt_ascension);
    *radius_vector = calcSunRadVector(T);
}

// Sun's topocentric (as seen from the observer's place on the Earth's surface) horizontal coordinates, in degrees
void calcHorizontalCoordinates(struct JulianDay jd, double latitude, double longitude,
                               double *azimuth, double *elevation)
{
    double T = calcJulianCent(jd);
    double GMST = calcGrMeanSiderealTime(jd);

    double ra, dec;
    calcSolarCoordinates(T, &ra, &dec);

    double H = GMST + longitude - ra;
    equatorial2horizontal(H, dec, latitude, azimuth, elevation);

    *azimuth += 180;  // measured from the North
    *elevation += calcRefraction(*elevation);
}

// Helper function
void calcRiseSetTimes(double m[3], struct JulianDay jd, double latitude, double longitude, double h0)
{
    double T = calcJulianCent(jd);
    double GMST = calcGrMeanSiderealTime(jd);

    double ra, dec;
    calcSolarCoordinates(T, &ra, &dec);

    // Local hour angle at sunrise or sunset (??NaN if body is circumpolar)
    double H0 = degrees(acos((sin(radians(h0)) - sin(radians(latitude)) * sin(radians(dec))) /
                        (cos(radians(latitude)) * cos(radians(dec)))));

    m[0] = jd.m + wrapTo180(ra - longitude - GMST) / 360;
    m[1] = m[0] - H0 / 360;
    m[2] = m[0] + H0 / 360;
}

// Find the times of sunrise, transit, and sunset, in hours
void calcSunriseSunset(struct JulianDay jd, double latitude, double longitude,
                       double *transit, double *sunrise, double *sunset, double altitude, int iterations)
{
    double m[3], times[3];
    m[0] = 0.5 - longitude / 360;

    for (int i = 0; i <= iterations; i++)
        for (int event = 0; event < 3; event++)
        {
            jd.m = m[event];
            calcRiseSetTimes(times, jd, latitude, longitude, altitude);
            m[event] = times[event];

            // First iteration, approximate rise and set times
            if (i == 0)
            {
                m[1] = times[1];
                m[2] = times[2];
                break;
            }
        }

    *transit = m[0] * 24;
    *sunrise = m[1] * 24;
    *sunset = m[2] * 24;
}

//======================================================================================================================
// Wrapper functions
//
// All calculations assume time inputs in Coordinated Universal Time (UTC)
//======================================================================================================================

void calcEquationOfTimeFromUnix(unsigned long utc, double *E)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcEquationOfTime(jd, E);
}

void calcEquationOfTimeFromDate(int year, int month, int day, int hour, int minute, int second, double *E)
{
    struct JulianDay jd = jdFromDate(year, month, day, hour, minute, second);
    calcEquationOfTime(jd, E);
}

void calcEquatorialCoordinatesFromUnix(unsigned long utc,
                                       double *rt_ascension, double *declination, double *radius_vector)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcEquatorialCoordinates(jd, rt_ascension, declination, radius_vector);
}

void calcEquatorialCoordinatesFromDate(int year, int month, int day, int hour, int minute, int second,
                                       double *rt_ascension, double *declination, double *radius_vector)
{
    struct JulianDay jd = jdFromDate(year, month, day, hour, minute, second);
    calcEquatorialCoordinates(jd, rt_ascension, declination, radius_vector);
}

void calcHorizontalCoordinatesFromUnix(unsigned long utc, double latitude, double longitude,
                                       double *azimuth, double *elevation)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcHorizontalCoordinates(jd, latitude, longitude, azimuth, elevation);
}

void calcHorizontalCoordinatesFromDate(int year, int month, int day, int hour, int minute, int second,
                                       double latitude, double longitude, double *azimuth, double *elevation)
{
    struct JulianDay jd = jdFromDate(year, month, day, hour, minute, second);
    calcHorizontalCoordinates(jd, latitude, longitude, azimuth, elevation);
}

void calcSunriseSunsetFromUnix(unsigned long utc, double latitude, double longitude,
                               double *transit, double *sunrise, double *sunset)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcSunriseSunset(jd, latitude, longitude, transit, sunrise, sunset, SUNRISESET_STD_ALTITUDE, 1);
}

void calcSunriseSunsetFromDate(int year, int month, int day, double latitude, double longitude,
                               double *transit, double *sunrise, double *sunset)
{
    struct JulianDay jd = jdFromDate(year, month, day, 0, 0, 0);
    calcSunriseSunset(jd, latitude, longitude, transit, sunrise, sunset, SUNRISESET_STD_ALTITUDE, 1);
}

void calcCivilDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                               double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, CIVIL_DAWNDUSK_STD_ALTITUDE, 1);
}

void calcCivilDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                               double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromDate(year, month, day, 0, 0, 0);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, CIVIL_DAWNDUSK_STD_ALTITUDE, 1);
}

void calcNauticalDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                                  double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, NAUTICAL_DAWNDUSK_STD_ALTITUDE, 1);
}

void calcNauticalDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                                  double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromDate(year, month, day, 0, 0, 0);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, NAUTICAL_DAWNDUSK_STD_ALTITUDE, 1);
}

void calcAstronomicalDawnDuskFromUnix(unsigned long utc, double latitude, double longitude,
                                      double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromUnix(utc);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE, 1);
}

void calcAstronomicalDawnDuskFromDate(int year, int month, int day, double latitude, double longitude,
                                      double *transit, double *dawn, double *dusk)
{
    struct JulianDay jd = jdFromDate(year, month, day, 0, 0, 0);
    calcSunriseSunset(jd, latitude, longitude, transit, dawn, dusk, ASTRONOMICAL_DAWNDUSK_STD_ALTITUDE, 1);
}
