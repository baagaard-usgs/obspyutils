#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#

import numpy

def _copyCoefs(coefs):
    coefCopy = []
    for coef in coefs:
        coefCopy.append(coef.copy())
    return coefCopy
    
def denoise(stream, remove_bg=True, preevent_window=10.0, preevent_threshold_reduction=2.0, store_noise=False, wavelet="coif4"):
    try:
        import pywt
    except ImportError:
        raise ImportError("_denoise() requires PyWavelets (pywt) Python module.")
    
    for tr in stream:
        tr.dataOrig = tr.data.copy()
        coefs = pywt.wavedec(tr.data, wavelet, mode="constant")
        tr.coefsOrig = _copyCoefs(coefs)

        coefsNoise = []

        if remove_bg:
            for coef in coefs:
                numCoef = coef.shape[-1]
                std = numpy.std(coef)
                mean = numpy.mean(coef)
                kurt = numpy.sum((coef-mean)**4) / (numCoef*std**4) - 3
                threshold = (24.0 / (numCoef*(1.0-0.9)))**0.5
                print("Preprocessing kurt: %f, threshold: %f" % (kurt, threshold,))
                mask = numpy.abs(kurt) <= threshold
                if mask:
                    coefsNoise.append(coef.copy())
                    coef *= 0.0
                else:
                    coefsNoise.append(0.0*coef)

        if preevent_window is not None and preevent_window > 0.0:
            numPtsPre = preevent_window*tr.stats.sampling_rate
            nlevels = len(coefs)-1
            for i,coef in enumerate(coefs[1:]):
                level = nlevels - i
                numCoefTarget = int(numPtsPre / 2**level)
                coefPre = coef[:numCoefTarget]
                numCoefPre = coefPre.shape[-1]
                median = numpy.median(numpy.abs(coefPre))
                std = median / 0.6745
                threshold = std * (2.0*numpy.log(numCoefPre)) / preevent_threshold_reduction
                print("Postprocessing threshold level: %d, : %f" % (level, threshold,))
                mask = numpy.abs(coef) < threshold
                coefsNoise[1+i][mask] += coef[mask]
                coefsNoise[1+i][~mask] += threshold*numpy.sign(coef[~mask])
                coef[mask] = 0.0
                coef[~mask] -= threshold*numpy.sign(coef[~mask])

        # Signal to noise ratio
        cArray,cSlices = pywt.coeffs_to_array(coefs)
        cArrayN,cSlices = pywt.coeffs_to_array(coefsNoise)
        mask = numpy.abs(cArray) > 0.0
        rmsSignal = numpy.sqrt(numpy.mean(cArray[mask]**2))
        rmsNoise = numpy.sqrt(numpy.mean(cArrayN[mask]**2))
        tr.StoN = rmsSignal/rmsNoise
        
        tr.data = pywt.waverec(coefs, wavelet, mode="constant")
        if tr.data.shape[-1] > tr.dataOrig.shape[-1]:
            tr.data = tr.data[:-1]
        if store_noise:
            tr.dataNoise = pywt.waverec(coefsNoise, wavelet, mode="constant")
            if tr.dataNoise.shape[-1] > tr.dataOrig.shape[-1]:
                tr.dataNoise = tr.dataNoise[:-1]
    return
            
