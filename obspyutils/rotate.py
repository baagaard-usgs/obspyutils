# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy
import numpy

# ----------------------------------------------------------------------
def _newTrace(data, header, component):
    """
    Create new trace from data, old header, and new component name.
    """
    channel = header.channel
    header.update({'channel': channel[:-1]+component})
    tr = obspy.core.Trace(data, header)
    return tr


# ----------------------------------------------------------------------
def toENZ(inventory, stream):

    from subset import streamByStation
    streamSt = streamByStation(stream)

    from math import sin,cos,pi
    
    tracesR = []
    for stkey,streamO in streamSt.items():
        ncode, scode = tuple(stkey.split('.'))
        station = inventory.select(network=ncode, station=scode).networks[0].stations[0]

        tr0 = streamO.traces[0]
        shape = tr0.data.shape
        dataE = numpy.zeros(shape)
        dataN = numpy.zeros(shape)
        dataZ = numpy.zeros(shape)

        starttime = tr0.stats.starttime
        endtime = tr0.stats.endtime 
        
        for tr in streamO.traces:
            cinventory = inventory.select(network=tr.stats.network, station=tr.stats.station, channel=tr.stats.channel)
            channel = cinventory.networks[0].stations[0].channels[0]
            trC = tr.copy().trim(starttime=starttime, endtime=endtime, pad=True, fill_value=0)
            azR = channel.azimuth.real * pi/180.0
            dipR = channel.dip.real * pi/180.0
            dataE += sin(azR)*cos(dipR)*trC.data
            dataN += cos(azR)*cos(dipR)*trC.data
            dataZ += -sin(dipR)*trC.data

        header = tr0.stats.copy()
        trE = _newTrace(dataE, header, 'E')
        trN = _newTrace(dataN, header, 'N')
        trZ = _newTrace(dataZ, header, 'Z')

        tracesR += [trE, trN, trZ]

    streamR = obspy.core.Stream(traces=tracesR)
    return streamR


# End of file