# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import math
import obspy

#-----------------------------------------------------------------------
def addLocation(inventory, stream):
    """
    Add longitude, latitude, and elevation to trace stats.
    """
    for trace in stream.traces:
        ist = inventory.select(network=trace.stats.network, station=trace.stats.station)
        if len(ist.networks) == 0 or len(ist.networks[0].stations) == 0:
            raise ValueError("Could not find station '%s.%s' in inventory." % trace.stats.network, trace.stats.station)
        stmeta = ist.networks[0].stations[0]
        info = {'longitude': stmeta.longitude,
                'latitude': stmeta.latitude,
                'elevation': stmeta.elevation,
                }
        trace.stats.update(info)
    return


#-----------------------------------------------------------------------
def addAzimuthDist(stream, epicenter, projection):
    """
    Add azimith and distance to stream stats.
    Epicenter is given in (longitude, latitude) and projection is a pyproj projection.
    
    projection = pyproj.Proj(proj='utm', zone=params.utmZone, ellps='WGS84')
    """
    # Project epicenter
    epicenterXY = projection(epicenter[0], epicenter[1])

    for trace in stream.traces:

        stXY = projection(trace.stats.longitude, trace.stats.latitude)
        epidist = ((stXY[0]-epicenterXY[0])**2 + (stXY[1]-epicenterXY[1])**2)**0.5
        azimuth = math.atan2(stXY[0]-epicenterXY[0], stXY[1]-epicenterXY[1])/math.pi*180.0
        if azimuth > 0:
            backazimuth = azimuth - 180.0
        else:
            backazimuth = azimuth + 180.0

        info = {'azimuth': azimuth,
                'back_azimuth': 180-azimuth,
                'distance': epidist}
        trace.stats.update(info)
    return


#-----------------------------------------------------------------------
def missing(inventory, stream):

    missing = {}
    found = {}
    for network in inventory.networks:
        for station in network.stations:
            key = "%s.%s" % (network.code, station.code)
            missing[key] = station

    for trace in stream.traces:
        key = "%s.%s" % (trace.stats.network, trace.stats.station)
        if key in missing.keys():
            missing.pop(key)
        if not key in found.keys():
            found[key] = station

    import collections

    print "Missing:",len(missing)
    missingO = collections.OrderedDict(sorted(missing.items()))
    for key, station in missingO.items():
        print key

    print "\n"
    print "Found:",len(found)
    foundO = collections.OrderedDict(sorted(found.items()))
    for key, station in foundO.items():
        print key

