# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import math
import obspy

#-----------------------------------------------------------------------
def addLocation(inventory, stream, epicenter=None, projection=None):
    """

    proj = pyproj.Proj(proj='utm', zone=params.utmZone, ellps='WGS84')

    """

    
    # Create metadata dictionary
    metadata = {}
    for network in inventory.networks:
        for station in network.stations:
            key = "%s.%s" % (network.code, station.code)
            metadata[key] = station

    # Project epicenter
    epicenterXY = proj(epicenter[0], epicenter[1])

    # Set metadata in traces
    if epicenter:
        for trace in stream.traces:
            key = "%s.%s" % (trace.stats.network, trace.stats.station)
            stmeta = metadata[key]

            stXY = proj(stmeta.longitude, stmeta.latitude)
            epidist = ((stXY[0]-epicenterXY[0])**2 + (stXY[1]-epicenterXY[1])**2)**0.5
            azimuth = math.atan2(stXY[0]-epicenterXY[0], stXY[1]-epicenterXY[1])/math.pi*180.0
            if azimuth > 0:
                backazimuth = azimuth - 180.0
            else:
                backazimuth = azimuth + 180.0

            info = {'longitude': stmeta.longitude,
                    'latitude': stmeta.latitude,
                    'elevation': stmeta.elevation,
                    'azimuth': azimuth,
                    'back_azimuth': 180-azimuth,
                    'distance': epidist}
            trace.stats.update(info)
    else:
        for trace in stream.traces:
            key = "%s.%s" % (trace.stats.network, trace.stats.station)
            stmeta = metadata[key]

            info = {'longitude': stmeta.longitude,
                    'latitude': stmeta.latitude,
                    'elevation': stmeta.elevation,
                    }
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
    missingO = collections.OrderedDict(sorted(missing.items()))
    for key, station in missingO.items():
        print key
    print "Missing:",len(missing)

    foundO = collections.OrderedDict(sorted(found.items()))
    for key, station in foundO.items():
        print key
    print "Found:",len(found)

