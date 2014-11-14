# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy
import numpy

# ======================================================================
# GNSCusp
class GNSCusp(object):
  """
  Python script for reading in a GNS CUSP file.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, filename=None):
    """
    Constructor.
    """
    self._headerTotalLines = 16 + 4 + 6
    self._headerAlphaLines = 16
    self._dataNumLines = 0
    self._dataNumValues = 0

    self.filename = filename
    self.header = {}
    self.azimuth = {}
    self.acc = None
    self.vel = None
    self.disp = None

    # Read the file.
    if self.filename is None:
      raise ValueError("Name of GNSCusp file has not been set.")

    fin = open(self.filename, 'r')
    lines = fin.readlines()
    fin.close()
    for component in xrange(3):
      self._parseHeader(lines, component)
    del lines

    self._readData()

    relTime = self.header['start_time'] - self.header['origin_time']
    if relTime.days < 0:
        relTime = self.header['origin_time'] - self.header['start_time']
        startTime = -relTime.seconds - relTime.microseconds/1e+6
    else:
      startTime = relTime.seconds + relTime.microseconds/1e+6

    (nsamples, ncomps) = self.acc.shape
    dt = self.header['dt']
    self.time = numpy.linspace(startTime, startTime + (nsamples-1)*dt,
                               nsamples)
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _parseHeader(self, lines, component):
    """
    Read first header in file.
    """
    header = {}

    i = component*(self._headerTotalLines+self._dataNumValues*self._dataNumLines)

    # String header
    i += 1
    # Line 2
    fields = lines[i].split(); i += 1
    header['id'] = fields[1]

    # Line 3
    header['name'] = lines[i].rstrip()

    i = self._headerAlphaLines + \
        component*(self._headerTotalLines+self._dataNumValues*self._dataNumLines)

    # Integer header
    # line 1
    fields1 = lines[i].split(); i += 1
    header['origin_time'] = datetime.datetime(int(fields1[0]),
                                              int(fields1[1]),
                                              int(fields1[2]),
                                              int(fields1[3]),
                                              int(fields1[4]),
                                              int(fields1[5])/10,
                                              int(int(fields1[5])%10*100000))

    # line 2
    fields2 = lines[i].split(); i += 1

    # line 3
    fields3 = lines[i].split(); i += 1
    self.azimuth[component] = float(fields3[7])
    header['epicentral_dist'] = float(fields3[9])

    # line 4
    fields4 = lines[i].split(); i += 1
    header['nsamples'] = int(fields4[3])

    if int(fields1[8]) != 0:
      header['start_time'] = datetime.datetime(int(fields1[8]),
                                               int(fields1[9]),
                                               int(fields2[8]),
                                               int(fields2[9]),
                                               int(fields4[8]),
                                               int(fields4[9])/1000,
                                               int(int(fields4[9])%1000)*1000)
    else:
      header['start_time']= header['origin_time']

    # Floating point heading
    # line 1
    fields = lines[i].split(); i += 1
    header['toSI'] = float(fields[7])*1.0e-3

    # line 2
    fields = lines[i].split(); i += 1
    lat = -float(fields[0])
    lon = float(fields[1])
    header['station_lonlatelev'] = numpy.array([lon, lat, 0.0], dtype=numpy.float64)

    # line 3
    fields = lines[i].split(); i += 1
    header['duration'] = float(fields[0])
    header['dt'] = float(fields[5])

    # line 4
    i += 1

    # line 5
    i += 1

    # line 6
    i += 1

    #print header
    #print self.azimuth

    if component == 0:
      self._dataNumLines = header['nsamples'] / 10
      if header['nsamples'] % 10 > 0: self._dataNumLines += 1
      totalDataLines = len(lines) - 3*self._headerTotalLines
      self._dataNumValues = totalDataLines / (3*self._dataNumLines)
      self.header = header
    return


  def _readData(self):
    """
    Read data from file.
    """
    raw = numpy.genfromtxt(self.filename, delimiter=[8]*10)

    # component 0
    dataStart = self._headerTotalLines
    acc0 = self._extract(raw, dataStart, self._dataNumLines)
    dataStart += self._dataNumLines
    if self._dataNumValues == 3:
      vel0 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines
      disp0 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines

    # component 1
    dataStart = self._headerTotalLines + \
        1*(self._headerTotalLines + self._dataNumValues*self._dataNumLines)
    acc1 = self._extract(raw, dataStart, self._dataNumLines)
    dataStart += self._dataNumLines
    if self._dataNumValues == 3:
      vel1 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines
      disp1 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines

    # component 2
    dataStart = self._headerTotalLines + \
        2*(self._headerTotalLines + self._dataNumValues*self._dataNumLines)
    acc2 = self._extract(raw, dataStart, self._dataNumLines)
    dataStart += self._dataNumLines
    if self._dataNumValues == 3:
      vel2 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines
      disp2 = self._extract(raw, dataStart, self._dataNumLines)
      dataStart += self._dataNumLines

    # Assemble components
    self.acc = numpy.array( (acc0, acc1, acc2) ).transpose()
    self.acc *= self.header['toSI']
    if self._dataNumValues == 3:
      self.vel = numpy.array( (vel0, vel1, vel2) ).transpose()
      self.vel *= self.header['toSI']
      self.disp = numpy.array( (disp0, disp1, disp2) ).transpose()
      self.disp *= self.header['toSI']

    # Rotate to E,N,U
    R = numpy.zeros((3, 3), dtype=numpy.float64)
    for i,az in enumerate(self.azimuth.values()):
      if az <= 360.0:
        azR = az / 180.0 * math.pi
        R[i,:] = (math.sin(azR), math.cos(azR), 0.0)
      else:
        R[i,:] = (0.0, 0.0, 1.0)
    self.acc = numpy.dot(self.acc, R)
    if not self.vel is None:
      self.vel = numpy.dot(self.vel, R)
    if not self.disp is None:
      self.disp = numpy.dot(self.disp, R)

    return


  def _extract(self, raw, start, size):
    values = raw[start:start+size,:]
    okay = values < 999999.8    
    return values[okay]


#-----------------------------------------------------------------------
def tostream(sites, filepattern, channelCode="HN"):
    """
    Convert GNS Cusp files to obspy stream.
    """

    tracesAcc = []
    tracesDisp = []
    tracesVel = []
    for site in sites:
        filename = filepattern % site
        cusp = GNSCusp(filename)

        for ic,component in enumerate(["E","N","Z"]):
            channel = "%s%s" % (channelCode, component)

            metadata = {'network': network,
                        'station': station,
                        'channel': channel,
                        'longitude': points[ipt,0],
                        'latitude': points[ipt,1],
                        'starttime': originTime+t[0],
                        'delta': dt,
                    }

            trace = obspy.core.Trace(data=cusp.acc[:,ic], header=metadata)
            tracesAcc.append(trace)
                
            trace = obspy.core.Trace(data=cusp.vel[:,ic], header=metadata)
            tracesVel.append(trace)
                
            trace = obspy.core.Trace(data=cusp.disp[:,ic], header=metadata)
            tracesDisp.append(trace)
                
    streamAcc = obspy.core.Stream(traces=tracesAcc)
    streamVel = obspy.core.Stream(traces=tracesVel)
    streamDisp = obspy.core.Stream(traces=tracesDisp)
    return streamAcc, streamVel, streamDisp


# End of file
