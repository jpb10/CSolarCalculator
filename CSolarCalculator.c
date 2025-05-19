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

#ifdef ARDUINO
#include <Arduino.h>
#else
#include <math.h>
#endif

#include "CSolarCalculator.h"

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

JulianDay jdFromUnix(unsigned long unix_ts)
{
    JulianDay jd;
    jd.JD = (unsigned long)(unix_ts / 86400) + 2440587.5;
    jd.m = (unix_ts % 86400) / 86400.0;
    return jd;
}

JulianDay jdFromDate(int year, int month, int day,
                     int hour, int minute, int second)
{
    JulianDay jd;
    /* Valid from 1901 to 2099 (Van Flandern & Pulkkinen, 1979) */
    jd.JD = 367.0 * year - (int)(7 * (year + (month + 9) / 12) / 4) +
            (int)(275 * month / 9) + day + 1721013.5;
    jd.m = (hour + minute / 60.0 + second / 3600.0) / 24.0;
    return jd;
}

/*
 * Intermediate calculations
 *
 * NB. Time T is measured in Julian centuries
 *     (36525 ephemeris days from the epoch J2000.0)
 */

#ifndef ARDUINO
double radians(double deg)
{
    return deg * 0.017453292519943296;
}

double degrees(double rad)
{
    return rad * 57.295779513082321;
}
#endif

double wrapTo360(double angle)
{
    angle = fmod(angle, 360);
    if (angle < 0) angle += 360;
    return angle; /* [0, 360) */
}

double wrapTo180(double angle)
{
    angle = wrapTo360(angle + 180);
    return angle - 180; /* [-180, 180) */
}

double calcJulianCent(JulianDay jd)
{
    return (jd.JD - 2451545 + jd.m) / 36525;
}

double calcGeomMeanLongSun(double T)
{
    return wrapTo360(280.46646 + T * 36000.76983); /* in degrees */
}

double calcGeomMeanAnomalySun(double T)
{
    return wrapTo360(357.52911 + T * 35999.05029); /* in degrees */
}

double calcSunEqOfCenter(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    return sin(radians(M)) * (1.914602 - 0.004817 * T) +
           sin(2 * radians(M)) * 0.019993; /* in degrees */
}

double calcSunRadVector(double T)
{
    double M = calcGeomMeanAnomalySun(T);
    return 1.00014 - 0.01671 * cos(radians(M)) -
           0.00014 * cos(2 * radians(M)); /* in AUs */
}

double calcMeanObliquityOfEcliptic(double T)
{
    return 23.4392911 - T * 0.0130042; /* in degrees */
}

/* Mean geocentric equatorial coordinates, accurate to ~1 arcminute */
void calcSolarCoordinates(double T,
                          double *ra, double *dec)
{
    double L0 = calcGeomMeanLongSun(T);
    double C = calcSunEqOfCenter(T);
    double L = L0 + C - 0.00569; /* corrected for aberration */

    double eps = calcMeanObliquityOfEcliptic(T);
    *ra = degrees(atan2(cos(radians(eps)) * sin(radians(L)), cos(radians(L))));
    *dec = degrees(asin(sin(radians(eps)) * sin(radians(L))));
}

double calcGrMeanSiderealTime(JulianDay jd)
{
    double GMST = wrapTo360(100.46061837 + 0.98564736629 * (jd.JD - 2451545));
    return wrapTo360(GMST + 360.985647 * jd.m); /* in degrees */
}

void equatorial2horizontal(double H, double dec, double lat,
                           double *az, double *el)
{
    double xhor = cos(radians(H)) * cos(radians(dec)) * sin(radians(lat)) -
                  sin(radians(dec)) * cos(radians(lat));
    double yhor = sin(radians(H)) * cos(radians(dec));
    double zhor = cos(radians(H)) * cos(radians(dec)) * cos(radians(lat)) +
                  sin(radians(dec)) * sin(radians(lat));

    *az = degrees(atan2(yhor, xhor));
    *el = degrees(atan2(zhor, sqrt(xhor * xhor + yhor * yhor)));
}

/* Hour angle at sunrise or sunset, returns NaN if circumpolar */
double calcHourAngleRiseSet(double dec, double lat, double h0)
{
    return degrees(acos((sin(radians(h0)) -
                         sin(radians(lat)) * sin(radians(dec))) /
                        (cos(radians(lat)) * cos(radians(dec)))));
}

/* Approximate atmospheric refraction correction, in degrees */
double calcRefraction(double elev)
{
    if (elev < -0.575)
        /* (Zimmerman, 1981) */
        return -20.774 / tan(radians(elev)) / 3600;
    else
        /* (Sæmundsson, 1986) */
        return 1.02 / tan(radians(elev + 10.3 / (elev + 5.11))) / 60;
}

/*
 * Solar calculator
 */

/* Equation of time, in minutes of time */
void calcEquationOfTime(JulianDay jd,
                        double *E)
{
    double T = calcJulianCent(jd);
    double L0 = calcGeomMeanLongSun(T);

    double ra, dec;
    calcSolarCoordinates(T, &ra, &dec);

    *E = 4 * wrapTo180(L0 - 0.00569 - ra);
}

/* Sun's (geocentric) equatorial coordinates, in degrees and AUs */
void calcEquatorialCoordinates(JulianDay jd,
                               double *rt_ascension, double *declination,
                               double *radius_vector)
{
    double T = calcJulianCent(jd);
    calcSolarCoordinates(T, rt_ascension, declination);

    *rt_ascension = wrapTo360(*rt_ascension);
    *radius_vector = calcSunRadVector(T);
}

/* Sun's (topocentric) horizontal coordinates, in degrees */
void calcHorizontalCoordinates(JulianDay jd, double latitude, double longitude,
                               double *azimuth, double *elevation)
{
    double T = calcJulianCent(jd);
    double GMST = calcGrMeanSiderealTime(jd);

    double ra, dec;
    calcSolarCoordinates(T, &ra, &dec);

    double H = GMST + longitude - ra;
    equatorial2horizontal(H, dec, latitude, azimuth, elevation);

    *azimuth += 180; /* measured from the North */
    *elevation += calcRefraction(*elevation);
}

/* Find the times of sunrise, transit, and sunset, in hours */
void calcSunriseSunset(JulianDay jd, double latitude, double longitude,
                       double *transit, double *sunrise, double *sunset,
                       double altitude, int iterations)
{
    double m[3];
    m[0] = 0.5 - longitude / 360;

    int i, event;
    for (i = 0; i <= iterations; i++)
        for (event = 0; event < 3; event++)
        {
            jd.m = m[event];
            double T = calcJulianCent(jd);
            double GMST = calcGrMeanSiderealTime(jd);

            double ra, dec;
            calcSolarCoordinates(T, &ra, &dec);

            double m0 = jd.m + wrapTo180(ra - longitude - GMST) / 360;
            double d0 = calcHourAngleRiseSet(dec, latitude, altitude) / 360;

            if (event == 0) m[0] = m0;
            if (event == 1 || i == 0) m[1] = m0 - d0;
            if (event == 2 || i == 0) m[2] = m0 + d0;
            if (i == 0) break;
        }

    *transit = m[0] * 24;
    *sunrise = m[1] * 24;
    *sunset = m[2] * 24;
}
