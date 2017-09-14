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
        coefs = pywt.wavedec(tr.data, wavelet)
        tr.coefsOrig = _copyCoefs(coefs)

        if store_noise:
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
                    if store_noise:
                        coefsNoise.append(coef.copy())
                    coef *= 0.0
                elif store_noise:
                    coefsNoise.append(0.0*coef)

        if preevent_window is not None and preevent_window > 0.0:
            numPtsPre = preevent_window*tr.stats.sampling_rate
            nlevels = len(coefs)-1
            for i,coef in enumerate(coefs[1:]):
                level = nlevels - i
                numCoefTarget = numPtsPre / 2**level
                coefPre = coef[:numCoefTarget]
                numCoefPre = coefPre.shape[-1]
                median = numpy.median(numpy.abs(coefPre))
                std = median / 0.6745
                threshold = std * (2.0*numpy.log(numCoefPre)) / preevent_threshold_reduction
                print("Postprocessing threshold level: %d, : %f" % (level, threshold,))
                mask = numpy.abs(coef) < threshold
                if store_noise:
                    coefsNoise[1+i][mask] += coef[mask]
                    coefsNoise[1+i][~mask] += threshold*numpy.sign(coef[~mask])
                coef[mask] = 0.0
                coef[~mask] -= threshold*numpy.sign(coef[~mask])

        tr.coefs = coefs
        tr.data = pywt.waverec(coefs, wavelet)
        if store_noise:
            tr.coefsNoise = coefsNoise
            tr.dataNoise = pywt.waverec(coefsNoise, wavelet)
    return
            
