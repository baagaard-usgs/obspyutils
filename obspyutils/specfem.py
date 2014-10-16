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
                        'longitude': longitude,
                        'latitude': latitude,
                        'starttime': originDateTime+t[0],
                        'delta': dt,
                        }

            trace = obspy.core.Trace(data=data, header=metadata)
            traces.append(trace)
                
    stream = obspy.core.Stream(traces=traces)
    return stream


# ----------------------------------------------------------------------
# ToObsPyApp
class ToObspyApp(object):
    """
    Python script for converting ASCII SPECFEM3D waveform output to ObsPy stream.
    """

    def __init__(self):
        """
        Constructor.
        """
        self.filenameIn = None
        self.filenameOut = None

        self.originTime = None
        self.epicenter = None
        self.utmZone = None

        self.channelCode = None
        self.dataType = None
        self.dataDir = None
        return


    def run(self):
        """
        Run conversion application.
        """
        import obspyutils.utils.io as io
        import obspyutils.utils.metadata as metadata

        s = tostream(self.filenameIn, self.dataDir, self.originTime, self.channelCode, self.dataType)

        # Add azimuth and distance
        if self.epicenter and self.utmZone:
            import pyproj
            projection = pyproj.Proj(proj='utm', zone=self.utmZone, ellps='WGS84')
            metadata.addAzimuthDist(s, self.epicenter, projection)
            
        io.pickle(self.filenameOut, s)

        return
  
#-----------------------------------------------------------------------
def writeCMT(event, centroidId="smi:nc.anss.org/origin/TMTS", hdur=0.0, filename="DATA/CMTSOLUTION", extractDC=False):
    """
    Write SPECFEM3D CMT file given event.
    """
    import event as eventutils
    import momenttensor

    evtname = eventutils.event_name(event)

    centroid = eventutils.find_origin(event, centroidId)
    time = centroid.time

    mt = eventutils.moment_tensor(event)
    Mw = momenttensor.Mw(mt)

    fout = open(filename, "w")
    fout.write("PDE %d %d %d %d %d %f" % (time.year, time.month, time.day, time.hour, time.minute, time.second))
    fout.write(" %.4f %.4f %.2f" % (centroid.latitude, centroid.longitude, centroid.depth/1000.0))
    fout.write(" %.2f %.2f" % (Mw, Mw))
    fout.write(" %s\n" % evtname)
    
    fout.write("event name: %s\n" % evtname)
    fout.write("time shift: %.1f\n" % 0.0)
    fout.write("half duration: %.1f\n" % hdur)
    fout.write("latitude: %.5f\n" % centroid.latitude)
    fout.write("longitude: %.5f\n" % centroid.longitude)
    fout.write("depth: %.3f\n" % (centroid.depth/1000.0))

    t = mt.tensor
    tscale = 1.0e+7
    if not extractDC:
        fout.write("Mrr: %12.4e\n" % (t.m_rr*tscale))
        fout.write("Mtt: %12.4e\n" % (t.m_tt*tscale))
        fout.write("Mpp: %12.4e\n" % (t.m_pp*tscale))
        fout.write("Mrt: %12.4e\n" % (t.m_rt*tscale))
        fout.write("Mrp: %12.4e\n" % (t.m_rp*tscale))
        fout.write("Mtp: %12.4e\n" % (t.m_tp*tscale))
    else:
        mDC = momenttensor.extractDC(mt)
        fout.write("Mrr: %12.4e\n" % (mDC[0,0]*tscale))
        fout.write("Mtt: %12.4e\n" % (mDC[0,1]*tscale))
        fout.write("Mpp: %12.4e\n" % (mDC[0,2]*tscale))
        fout.write("Mrt: %12.4e\n" % (mDC[1,1]*tscale))
        fout.write("Mrp: %12.4e\n" % (mDC[1,2]*tscale))
        fout.write("Mtp: %12.4e\n" % (mDC[2,2]*tscale))

    return

# End of file
