# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy
import numpy

#-----------------------------------------------------------------------
def tostream(filename, originTime=None, projection=None, channelCode="HX", dataType='vel'):
    """
    Convert HDF5 waveform output from PyLith simulation to obspy stream.
    """
    import h5py

    h5 = h5py.File(filename, "r", driver='sec2')
    points = h5['/geometry/vertices'][:]
    if dataType == 'disp':
        field = 'displacement'
    elif dataType == 'vel':
        field = 'velocity'
    elif dataType == 'acc':
        field = 'acceleration'
    else:
        raise ValueError("Unknown data type '%s'." % dataType)
    t = h5['/time']
    stations = h5['stations']
    data = h5['/vertex_fields/%s' % field][:]
    h5.close()

    if projection:
        lon,lat = projection(points[:,0], points[:,1], inverse=True)
        points[:,0] = lon
        points[:,1] = lat

    npts, ndims = points.shape
    
    nsteps, npts2, ncomps = data.shape
    assert(npts == npts2)
    assert(3 == ndims)
    assert(3 == ncomps)

    dt = t[1]-t[0] # Assume uniform time step

    for ipt in xrange(npts):

        (network, station) = stations[ipt].split('.')

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

            trace = obspy.core.Trace(data=data[:,ipt,ic], header=metadata)
            traces.append(trace)
                
    stream = obspy.core.Stream(traces=traces)
    return stream


# End of file
