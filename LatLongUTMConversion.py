#!/usr/bin/env python
# Lat Long - UTM, UTM - Lat Long conversions

# original source: http://webapp.decouto.bm/media/bng/about-bng.html
# minor changes were made by @muloka, mainly linting this source

from math import pi, sin, cos, tan, sqrt

#LatLong- UTM conversion..h
#definitions for lat/long to UTM and UTM to lat/lng conversions
#include <string.h>

_deg2rad = pi / 180.0
_rad2deg = 180.0 / pi

_EquatorialRadius = 2
_eccentricitySquared = 3

_ellipsoid = [
# id, Ellipsoid name, Equatorial Radius, square of eccentricity
# first once is a placeholder only, To allow array indices to match id numbers
    [-1, "Placeholder", 0, 0],
    [1, "Airy", 6377563, 0.00667054],
    [2, "Australian National", 6378160, 0.006694542],
    [3, "Bessel 1841", 6377397, 0.006674372],
    [4, "Bessel 1841 (Nambia] ", 6377484, 0.006674372],
    [5, "Clarke 1866", 6378206, 0.006768658],
    [6, "Clarke 1880", 6378249, 0.006803511],
    [7, "Everest", 6377276, 0.006637847],
    [8, "Fischer 1960 (Mercury] ", 6378166, 0.006693422],
    [9, "Fischer 1968", 6378150, 0.006693422],
    [10, "GRS 1967", 6378160, 0.006694605],
    [11, "GRS 1980", 6378137, 0.00669438],
    [12, "Helmert 1906", 6378200, 0.006693422],
    [13, "Hough", 6378270, 0.00672267],
    [14, "International", 6378388, 0.00672267],
    [15, "Krassovsky", 6378245, 0.006693422],
    [16, "Modified Airy", 6377340, 0.00667054],
    [17, "Modified Everest", 6377304, 0.006637847],
    [18, "Modified Fischer 1960", 6378155, 0.006693422],
    [19, "South American 1969", 6378160, 0.006694542],
    [20, "WGS 60", 6378165, 0.006693422],
    [21, "WGS 66", 6378145, 0.006694542],
    [22, "WGS-72", 6378135, 0.006694318],
    [23, "WGS-84", 6378137, 0.00669438]
]

#Reference ellipsoids derived from Peter H. Dana's website-
#http://www.utexas.edu/depts/grg/gcraft/notes/datum/elist.html
#Department of Geography, University of Texas at Austin
#Internet: pdana@mail.utexas.edu
#3/22/95

#Source
#Defense Mapping Agency. 1987b. DMA Technical Report: Supplement to Department of Defense World Geodetic System
#1984 Technical Report. Part I and II. Washington, DC: Defense Mapping Agency

#def LLtoUTM(int ReferenceEllipsoid, const double Lat, const double Long,
#			 double &UTMNorthing, double &UTMEasting, char* UTMZone)


def LLtoUTM(ReferenceEllipsoid, Lat, Long, zone=None):
    """converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees
    Written by Chuck Gantz- chuck.gantz@globalstar.com"""

    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    k0 = 0.9996

    #Make sure the longitude is between -180.00 .. 179.9
    LongTemp = (Long + 180) - int((Long + 180) / 360) * 360 - 180  # -180.00 .. 179.9

    LatRad = Lat * _deg2rad
    LongRad = LongTemp * _deg2rad

    if zone is None:
        ZoneNumber = int((LongTemp + 180) / 6) + 1
    else:
        ZoneNumber = zone

    if Lat >= 56.0 and Lat < 64.0 and LongTemp >= 3.0 and LongTemp < 12.0:
        ZoneNumber = 32

    # Special zones for Svalbard
    if Lat >= 72.0 and Lat < 84.0:
        if LongTemp >= 0.0 and LongTemp < 9.0:
            ZoneNumber = 31
        elif LongTemp >= 9.0 and LongTemp < 21.0:
            ZoneNumber = 33
        elif LongTemp >= 21.0 and LongTemp < 33.0:
            ZoneNumber = 35
        elif LongTemp >= 33.0 and LongTemp < 42.0:
            ZoneNumber = 37

    LongOrigin = (ZoneNumber - 1) * 6 - 180 + 3  # +3 puts origin in middle of zone
    LongOriginRad = LongOrigin * _deg2rad

    #compute the UTM Zone from the latitude and longitude
    UTMZone = "%d%c" % (ZoneNumber, _UTMLetterDesignator(Lat))

    eccPrimeSquared = (eccSquared) / (1 - eccSquared)
    N = a / sqrt(1 - eccSquared * sin(LatRad) * sin(LatRad))
    T = tan(LatRad) * tan(LatRad)
    C = eccPrimeSquared * cos(LatRad) * cos(LatRad)
    A = cos(LatRad) * (LongRad - LongOriginRad)

    M = a * ((1
            - eccSquared / 4
            - 3 * eccSquared * eccSquared / 64
            - 5 * eccSquared * eccSquared * eccSquared / 256) * LatRad
           - (3 * eccSquared / 8
              + 3 * eccSquared * eccSquared / 32
              + 45 * eccSquared * eccSquared * eccSquared / 1024) * sin(2 * LatRad)
           + (15 * eccSquared * eccSquared / 256 + 45 * eccSquared * eccSquared * eccSquared / 1024) * sin(4 * LatRad)
           - (35 * eccSquared * eccSquared * eccSquared / 3072) * sin(6 * LatRad))

    UTMEasting = (k0 * N * (A + (1 - T + C) * A * A * A / 6
                        + (5 - 18 * T + T * T + 72 * C - 58 * eccPrimeSquared) * A * A * A * A * A / 120)
                  + 500000.0)

    UTMNorthing = (k0 * (M + N * tan(LatRad) * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A * A * A * A / 24
                                        + (61
                                           - 58 * T
                                           + T * T
                                           + 600 * C
                                           - 330 * eccPrimeSquared) * A * A * A * A * A * A / 720)))

    if Lat < 0:
        UTMNorthing = UTMNorthing + 10000000.0  # 10000000 meter offset for southern hemisphere
    return (UTMZone, UTMEasting, UTMNorthing)


def _UTMLetterDesignator(Lat):
    """This routine determines the correct UTM letter designator for the given latitude
    returns 'Z' if latitude is outside the UTM limits of 84N to 80S
    Written by Chuck Gantz- chuck.gantz@globalstar.com"""

    if 84 >= Lat >= 72:
        return 'X'
    elif 72 > Lat >= 64:
        return 'W'
    elif 64 > Lat >= 56:
        return 'V'
    elif 56 > Lat >= 48:
        return 'U'
    elif 48 > Lat >= 40:
        return 'T'
    elif 40 > Lat >= 32:
        return 'S'
    elif 32 > Lat >= 24:
        return 'R'
    elif 24 > Lat >= 16:
        return 'Q'
    elif 16 > Lat >= 8:
        return 'P'
    elif 8 > Lat >= 0:
        return 'N'
    elif 0 > Lat >= -8:
        return 'M'
    elif -8 > Lat >= -16:
        return 'L'
    elif -16 > Lat >= -24:
        return 'K'
    elif -24 > Lat >= -32:
        return 'J'
    elif -32 > Lat >= -40:
        return 'H'
    elif -40 > Lat >= -48:
        return 'G'
    elif -48 > Lat >= -56:
        return 'F'
    elif -56 > Lat >= -64:
        return 'E'
    elif -64 > Lat >= -72:
        return 'D'
    elif -72 > Lat >= -80:
        return 'C'
    else:
        return 'Z'  # if the Latitude is outside the UTM limits

#void UTMtoLL(int ReferenceEllipsoid, const double UTMNorthing, const double UTMEasting, const char* UTMZone,
#			  double& Lat,  double& Long )


def UTMtoLL(ReferenceEllipsoid, northing, easting, zone,
            false_E=None, false_N=None, orig_lat=None, orig_lon=None,
            k0_scale=None):
    """converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
    East Longitudes are positive, West longitudes are negative.
    North latitudes are positive, South latitudes are negative
    Lat and Long are in decimal degrees.
    Written by Chuck Gantz- chuck.gantz@globalstar.com
    Converted to Python by Russ Nelson <nelson@crynwr.com>"""

    k0 = 0.9996
    a = _ellipsoid[ReferenceEllipsoid][_EquatorialRadius]
    eccSquared = _ellipsoid[ReferenceEllipsoid][_eccentricitySquared]
    e1 = (1 - sqrt(1 - eccSquared)) / (1 + sqrt(1 - eccSquared))
    #NorthernHemisphere; //1 for northern hemispher, 0 for southern

    if false_E is None:
        false_E = 500000.0  # standard UTM false easting is 500,000 meters
    x = easting - false_E  # remove false easting offset for longitude
    y = northing

    if false_N is None:
        ZoneLetter = zone[-1]
        ZoneNumber = int(zone[:-1])
        if ZoneLetter <= 'N':
            y -= 10000000.0         # remove 10,000,000 meter offset used for southern hemisphere
            # NorthernHemisphere = 0  # point is in southern hemisphere
        # else:
            # NorthernHemisphere = 1  # point is in northern hemisphere

        LongOrigin = (ZoneNumber - 1) * 6 - 180 + 3  # +3 puts origin in middle of zone
    else:
        # Assume other UTM params are available
        # if orig_lat >= 0:
        #    NorthernHemisphere = 1
        # else:
        #    NorthernHemisphere = 0
        y -= false_N
        LongOrigin = orig_lon
        k0 = k0_scale

    eccPrimeSquared = (eccSquared) / (1 - eccSquared)
    e2 = eccSquared

    if orig_lat is None:
        M0 = 0
    else:
        lat0Rad = orig_lat * _deg2rad

        M0 = a * (
            (1
                 - e2 / 4
                 - 3 * e2 * e2 / 64
                 - 5 * e2 * e2 * e2 / 256) * lat0Rad
            - (3 * e2 / 8
               + 3 * e2 * e2 / 32
               + 45 * e2 * e2 * e2 / 1024) * sin(2 * lat0Rad)
            + (15 * e2 * e2 / 256
               + 45 * e2 * e2 * e2 / 1024) * sin(4 * lat0Rad)
            - (35 * e2 * e2 * e2 / 3072) * sin(6 * lat0Rad)
            )

    M = M0 + y / k0
    mu = M / (a * (1 - eccSquared / 4 - 3 * eccSquared * eccSquared / 64 - 5 * eccSquared * eccSquared * eccSquared / 256))

    phi1Rad = (mu + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * sin(2 * mu)
               + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * sin(4 * mu)
               + (151 * e1 * e1 * e1 / 96) * sin(6 * mu)
               + (1097 * e1 * e1 * e1 * e1 / 512) * sin(8 * mu)  # XXX ???
               )

    # phi1 = phi1Rad * _rad2deg;

    N1 = a / sqrt(1 - eccSquared * sin(phi1Rad) * sin(phi1Rad))
    T1 = tan(phi1Rad) * tan(phi1Rad)
    C1 = eccPrimeSquared * cos(phi1Rad) * cos(phi1Rad)
    R1 = a * (1 - eccSquared) / pow(1 - eccSquared * sin(phi1Rad) * sin(phi1Rad), 1.5)
    D = x / (N1 * k0)

    Lat = phi1Rad - (N1 * tan(phi1Rad) / R1) * \
                        (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eccPrimeSquared) * D * D * D * D / 24
                       + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * eccPrimeSquared - 3 * C1 * C1) * D * D * D * D * D * D / 720)
    Lat = Lat * _rad2deg

    Long = (D - (1 + 2 * T1 + C1) * D * D * D / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * eccPrimeSquared + 24 * T1 * T1)
            * D * D * D * D * D / 120) / cos(phi1Rad)
    Long = LongOrigin + Long * _rad2deg
    return (Lat, Long)


def bng_to_latlon(east, north):
    false_E = 550000.0
    false_N = 100000.0

    origin_lat = 32.0
    origin_lon = -64.75
    scale = 1.0

    ref_ellipsoid = 23  # WGS84
    lat, lon = UTMtoLL(ref_ellipsoid, north, east, None,
                       false_E=false_E, false_N=false_N, k0_scale=scale,
                       orig_lat=origin_lat, orig_lon=origin_lon)
    return lat, lon


def main():

    # Coords of High Point house
    hp_lat = 32.249546
    hp_lon = -64.855876
    hp_N = 127672
    hp_E = 540046

    false_E = 550000.0
    false_N = 100000.0

    origin_lat = 32.0
    origin_lon = -64.75
    scale = 1.0

    for ref_ellipsoid in [23]:  # WGS84
        lat, lon = UTMtoLL(ref_ellipsoid, hp_N, hp_E, None,
                           false_E=false_E, false_N=false_N, k0_scale=scale,
                           orig_lat=origin_lat, orig_lon=origin_lon)

        print
        print "Ellipsoid", ref_ellipsoid
        print "Expect %8.5f, %9.5f" % (hp_lat, hp_lon)
        print "Got    %8.5f, %9.5f" % (lat, lon)

    print
    print
    print "Spreadsheet test"
    ss_N = 3569608.462
    ss_E = 325165.7693
    lat, lon = UTMtoLL(ref_ellipsoid, ss_N, ss_E, '20S')
    print "Expect %8.5f, %9.5f" % (hp_lat, hp_lon)
    print "Got    %8.5f, %9.5f" % (lat, lon)


if __name__ == '__main__':
    main()
