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

            metadata = {'network': network,
                        'station': station,
                        'channel': channel,
                        'longitude': longitude,
                        'latitude': latitude,
                        'starttime': originTime+t[0],
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
        import obspyutils.pickle as pickle
        import obspyutils.metadata as metadata

        s = tostream(self.filenameIn, self.dataDir, self.originTime, self.channelCode, self.dataType)

        # Add azimuth and distance
        if self.epicenter and self.utmZone:
            import pyproj
            projection = pyproj.Proj(proj='utm', zone=self.utmZone, ellps='WGS84')
            metadata.addAzimuthDist(s, self.epicenter, projection)
            
        pickle.pickle(self.filenameOut, s)

        return
  
#-----------------------------------------------------------------------
def writeCMT(event, originId=None, mechanismId=None, hdur=0.0, filename="DATA/CMTSOLUTION"):
    """
    Write SPECFEM3D CMT file given event.
    """
    import event as eventutils
    import momenttensor

    evtname = eventutils.event_name(event)

    if originId is None:
        originId = "smi:nc.anss.org/origin/HYP2000"
    origin = eventutils.find_origin(event, originId)
    time = origin.time

    if mechanismId is None:
        mechanismId = "smi:nc.anss.org/momentTensor/TMTS"
    mt = eventutils.find_momenttensor(event, mechanismId)
    Mw = momenttensor.Mw(mt)

    fout = open(filename, "w")
    fout.write("PDE %5d%3d%3d%3d%3d%6.2f" % (time.year, time.month, time.day, time.hour, time.minute, time.second))
    fout.write(" %.4f %.4f %.2f" % (origin.latitude, origin.longitude, origin.depth/1000.0))
    fout.write(" %.2f %.2f" % (Mw, Mw))
    fout.write(" %s\n" % evtname)
    
    fout.write("event name: %s\n" % evtname)
    fout.write("time shift: %.1f\n" % 0.0)
    fout.write("half duration: %.1f\n" % hdur)
    fout.write("latitude: %.5f\n" % origin.latitude)
    fout.write("longitude: %.5f\n" % origin.longitude)
    fout.write("depth: %.3f\n" % (origin.depth/1000.0))

    t = mt.tensor
    tscale = 1.0e+7
    fout.write("Mrr: %12.4e\n" % (t.m_rr*tscale))
    fout.write("Mtt: %12.4e\n" % (t.m_tt*tscale))
    fout.write("Mpp: %12.4e\n" % (t.m_pp*tscale))
    fout.write("Mrt: %12.4e\n" % (t.m_rt*tscale))
    fout.write("Mrp: %12.4e\n" % (t.m_rp*tscale))
    fout.write("Mtp: %12.4e\n" % (t.m_tp*tscale))

    return


#-----------------------------------------------------------------------
def writeStations(inventory, filename="DATA/STATIONS"):
    """
    Write SPECFEM3D STATIONS file given station inventory.
    """
    depth = 0.0

    fout = open(filename, "w")
    for network in inventory.networks:
        for station in network.stations:
            fout.write("%s %s %.4f %.4f %.1f %.1f\n" %\
                       (station.code, network.code, station.latitude, station.longitude, station.elevation, depth))
    fout.close()
    return


# End of file
