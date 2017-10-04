# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import numpy
import obspy
import scipy.optimize
import scipy.interpolate
import scipy.integrate

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
    p = scipy.optimize.curve_fit(_linearvel, tfit, vfit)
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
def baseline_correction_spline(stream, window=2.0):
    for tr in stream.traces:
        t = tr.times()
        knotsN = int(0.1*window*tr.stats.sampling_rate)
        knots = numpy.linspace(t[knotsN/2], t[-knotsN/2], len(t)/knotsN)
        spline = scipy.interpolate.LSQUnivariateSpline(tr.times(), tr.data, t=knots, k=5)
        tr.data -= spline(t)
    return


# ----------------------------------------------------------------------
def correction_preonly(data, preevent_window=10.0):
        if isinstance(data, obspy.core.Stream):
            for tr in data:
                correction_preonly(tr, preevent_window)
        else:
            dt = data.stats.delta
            numPreWindow = 1+int(preevent_window/dt)

            poly = numpy.polyfit(data.times()[:numPreWindow], scipy.integrate.cumtrapz(data.data[:numPreWindow], dx=dt, initial=0.0), deg=1)
            data.data -= poly[0]
        return

# ----------------------------------------------------------------------
def correction_constant(data):
        if isinstance(data, obspy.core.Stream):
            for tr in data:
                correction_constant(tr)
        else:
            dt = data.stats.delta
            vel = scipy.integrate.cumtrapz(data.data, dx=dt, initial=0.0)
            disp = scipy.integrate.cumtrapz(vel, dx=dt, initial=0.0)
            poly = numpy.polyfit(data.times(), disp, deg=2)
            data.data -= 2.0*poly[0]
        return

# ----------------------------------------------------------------------
def integrate_acc(data):
        """
        """
        if isinstance(data, obspy.core.Stream):
            tracesVel = []
            tracesDisp = []
            for tr in stream:
                trVel,trDisp = integrate_acc(tr)
                tracesVel.append(trVel)
                tracesDisp.append(trDisp)
            return (obspy.core.Stream(traces=tracesVel), obspy.core.Stream(traces=tracesDisp),)
        else:
            dt = data.stats.delta
            t = data.times()
            vel = scipy.integrate.cumtrapz(data.data, dx=dt, initial=0.0)
            disp = scipy.integrate.cumtrapz(vel, dx=dt, initial=0.0)
            poly = numpy.polyfit(t, disp, deg=2)
            # poly[0] should be 0 if we have a good baseline correction
            
            vel -= 2.0*poly[0]*t + poly[1]
            disp -= poly[0]*t**2 + poly[1]*t + poly[2]

            trVel = data.copy()
            trVel.data = vel

            trDisp = data.copy()
            trDisp.data = disp

            return (trVel, trDisp,)
            
        return


# ----------------------------------------------------------------------
def integrate_vel(data):
        """
        """
        if isinstance(data, obspy.core.Stream):
            tracesDisp = []
            for tr in stream:
                trDisp = integrate_vel(tr)
                tracesDisp.append(trDisp)
            return obspy.core.Stream(traces=tracesDisp)
        else:
            dt = data.stats.delta
            t = data.times()
            disp = scipy.integrate.cumtrapz(data.data, dx=dt, initial=0.0)
            poly = numpy.polyfit(t, disp, deg=1)
            
            disp -= 2.0*poly[0]*t + poly[1]
            disp -= poly[0]

            trDisp = data.copy()
            trDisp.data = disp

            return trDisp
            
        return


# End of file
