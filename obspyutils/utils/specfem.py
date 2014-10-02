# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy
import numpy

#-----------------------------------------------------------------------
def tostream(filename="DATA/STATIONS_FILTERED", dataDir="OUTPUT_FILES", originTime=None, channelCode="HX", dataType='vel'):
    """
    Collect ASCII waveform output from SPECFEM3D simulation and convert them to obspy stream.
    """
    fin = open(filename, "r")
    lines = fin.readlines()
    fin.close()

    traces = []
    for line in lines:
        fields = line.split()
        if len(fields) != 6:
            raise IOError("Unrecognized format for stations file.\nLine: '%s'") % line
        
        station = fields[0].strip()
        network = fields[1].strip()
        latitude = float(fields[2])
        longitude = float(fields[3])
        elevation = float(fields[4])

        if dataType[0:4] == "disp":
            suffix = "semd"
        elif dataType[0:3] == "vel":
            suffix = "semv"
        elif dataType[0:3] == "acc":
            suffix = "sema"
            
        for component in ["E","N","Z"]:
            channel = "%s%s" % (channelCode, component)
            wfilename = "%s/%s.%s.%s.%s" % \
              (dataDir, station, network, channel, suffix)
            raw = numpy.loadtxt(wfilename)
            t = raw[:,0]
            data = raw[:,1]
            dt = t[1]-t[0]

            try:
                originDateTime = obspy.core.UTCDateTime(originTime)
            except TypeError:
                raise TypeError("Cannot parse originTime '%s' into a UTCDateTime object." % originTime)
            
            metadata = {'network': network,
                        'station': station,
                        'channel': channel,
                        'starttime': originDateTime,
                        'delta': dt,
                        }

            trace = obspy.core.Trace(data=data, header=metadata)
            traces.append(trace)
                
    stream = obspy.core.Stream(traces=traces)
    return stream


# End of file
