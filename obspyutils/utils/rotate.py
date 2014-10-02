# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

#-----------------------------------------------------------------------
def rotate(stream, method="NE->RT"):
    from subset import groupByStation
    stations = groupByStation(stream)

    tracesR = []
    for stationO in stations.values():
        try:
            stationR = stationO.rotate('NE->RT')
            tracesR += stationR.traces
        except ValueError:
            print "Bad station.",stationO
    streamR = obspy.core.Stream(traces=tracesR)
    return streamR


# End of file
