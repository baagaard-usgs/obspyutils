# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import numpy
import obspy
from scipy import optimize
from scipy import interpolate

# ----------------------------------------------------------------------
def _linearvel(t, v0, af):
    t2 = t[-1]-30.0
    return v0+af*(t-t2)


# ----------------------------------------------------------------------
def _correctionV0(acc, tpre, ttail):

    # Remove pre-event mean
    ipre = numpy.argmax(acc.times() >= tpre)
    acc.data -= numpy.mean(acc.data[0:ipre])
    
    vel = acc.copy()
    vel.integrate()
    t = vel.times() - tpre

    # Curve fit line to tail of velocity record
    t2 = t[-1] - ttail
    mask = t >= t2
    tfit = t[mask]
    vfit = vel.data[mask]
    p = optimize.curve_fit(_linearvel, tfit, vfit)
    (v0, af) = p[0]

    # Average acceleration in strong part of record
    i1 = numpy.argmax(numpy.abs(acc.data) >= 0.05)
    t1 = t[i1]
    am = v0 / (t2-t1)
    
    i2 = numpy.argmax(t >= t2)

    mask = numpy.bitwise_and(t > t1, t <= t2)
    acc.data[mask] -= am
    
    mask = t > t2
    acc.data[mask] -= af
    
    return
    

# ----------------------------------------------------------------------
def baseline_correction_v0(stream, originTime, ttail=None):
    for tr in stream.traces:
        tpre = originTime - tr.stats.starttime
        if not ttail:
            ttail = tpre
        _correctionV0(tr, tpre, ttail)
    return


# ----------------------------------------------------------------------
def baseline_correction_spline(stream, knotsN=400):
    for tr in stream.traces:
        t = tr.times()
        knots = numpy.linspace(t[knotsN/2], t[-knotsN/2], len(t)/knotsN)
        spline = interpolate.LSQUnivariateSpline(tr.times(), tr.data, t=knots, k=5)
        tr.data -= spline(t)
    return


# End of file
