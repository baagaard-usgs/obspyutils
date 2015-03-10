# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================


import numpy
from obspy.core.event import MomentTensor,Tensor

#-----------------------------------------------------------------------
def toTensor(a):
    t = Tensor(m_rr=a[0,0], m_tt=a[1,1], m_pp=a[2,2],
               m_rt=a[0,1], m_rp=a[0,2], m_tp=a[1,2])
    return t


#-----------------------------------------------------------------------
def Mw(mt):
    from math import log10
    Mo = mt.scalar_moment*1.0e+7 # Nm -> dyne-cm
    return log10(Mo)/1.5-10.7


#-----------------------------------------------------------------------
def seismicMoment(Mw):
    """
    Seismic moment in Nm.
    """
    return 10**(1.5*(Mw+10.7)-7)


#-----------------------------------------------------------------------
def extractDC(mt, rescale=False):
    """
    WARNING: NOT VERIFIED!!!
    """
    t = mt.tensor
    M = numpy.array([[t.m_rr, t.m_rt, t.m_rp],
                     [t.m_rt, t.m_tt, t.m_tp],
                     [t.m_rp, t.m_tp, t.m_pp]],
                    dtype=numpy.float64)

    w,v = numpy.linalg.eig(M)
    wv = sorted(zip(w,v), key=lambda v: v[0])

    s3,v3 = wv[2]
    s2,v2 = wv[1]
    s1,v1 = wv[0]

    R = numpy.vstack((v1,v2,v3)).transpose()

    MdcP = numpy.array([[0.5*(s1-s3), 0, 0], [0, 0, 0], [0, 0, -0.5*(s1-s3)]],
                       dtype=numpy.float64)

    Mdc = numpy.dot(numpy.dot(R, MdcP), R.transpose())
    if rescale:
        MoDC = 0.5**0.5*numpy.tensordot(Mdc,Mdc)**0.5
        Mdc *= mt.scalar_moment / MoDC
    mtdc = mt.copy()
    mtdc.update({'tensor': toTensor(Mdc), 'scalar_moment': 0.5**0.5*numpy.tensordot(Mdc,Mdc)**0.5})
    #return mtdc
    return # Not verified


#-----------------------------------------------------------------------
def rescale(mt, Mo=None):
    """
    Rescale moment tensor to match seismic moment.
    """

    t = mt.tensor
    M = numpy.array([[t.m_rr, t.m_rt, t.m_rp],
                     [t.m_rt, t.m_tt, t.m_tp],
                     [t.m_rp, t.m_tp, r.m_pp]],
                    dtype=numpy.float64)

    Mo = Mo or mt.scalar_moment
    M *= Mo / (0.5**0.5*numpy.tensordot(M,M)**0.5)
    t.m_rr = M[0,0]
    t.m_rt = M[0,1]
    t.m_rp = M[0,2]
    t.m_tt = M[1,1]
    t.m_tp = M[1,2]
    t.m_pp = M[2,2]
    return

#-----------------------------------------------------------------------
def anglesToMT(strike_d, dip_d, rake_d, Mo):
    """
    Convert strike, dip, rake, to MT in Up,South,East coordinate system.
    """
    from math import pi,sin,cos
    strike = strike_d / 180.0 * pi
    dip = dip_d / 180.0 * pi
    rake = rake_d / 180.0 * pi

    Mxx = -Mo*(sin(dip)*cos(rake)*sin(2*strike) + sin(2*dip)*sin(rake)*sin(strike)**2)
    Myy = +Mo*(sin(dip)*cos(rake)*sin(2*strike) - sin(2*dip)*sin(rake)*cos(strike)**2)
    Mzz = +Mo*(sin(2*dip)*sin(rake))

    Mxy = +Mo*(sin(dip)*cos(rake)*cos(2*strike) + 0.5*sin(2*dip)*sin(rake)*sin(2*strike))
    Myz = -Mo*(cos(dip)*cos(rake)*sin(strike) - cos(2*dip)*sin(rake)*cos(strike))
    Mxz = -Mo*(cos(dip)*cos(rake)*cos(strike) + cos(2.0*dip)*sin(rake)*sin(strike))

    M = numpy.array([[+Mzz, +Mxz, -Myz],
                     [+Mxz, +Mxx, -Mxy],
                     [-Myz, -Mxy, +Myy]],
                    dtype=numpy.float64)

    mt = MomentTensor(scalar_moment=Mo, tensor=toTensor(M))

    return mt



