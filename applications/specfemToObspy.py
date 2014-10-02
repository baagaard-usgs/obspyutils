#!/usr/bin/env python
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================
#
# Python script for converting ASCII SPECFEM3D waveform output to ObsPy stream.

# ======================================================================
# ConvertApp
class ConvertApp(object):
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
        self.channelCode = None
        self.dataType = None
        self.dataDir = None
        return


    def run(self):
        """
        Run conversion application.
        """
        import obspyutils.utils.specfem as specfem
        import obspyutils.utils.io as io

        s = specfem.tostream(self.filenameIn, self.dataDir, self.originTime, self.channelCode, self.dataType)
        io.pickle(self.filenameOut, s)

        return
  

# ======================================================================
if __name__ == '__main__':

    usage = "%prog [options]\n\n" \
        "Convert ASCII SPECFEM3D waveform output to ObsPy stream."

    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    parser.add_option("-i", "--filein", dest="filenameIn",
                      type="string", metavar="FILE",
                      help="Read station list from FILE. [DATA/STATIONS_FILTERED]",
                      default="DATA/STATIONS_FILTERED")
    parser.add_option("-o", "--fileout", dest="filenameOut",
                      type="string", metavar="FILE",
                      help="Pickle stream to FILE. [waveforms_syn.p]",
                      default="waveforms_syn.p")
    parser.add_option("-s", "--origintime", dest="originTime",
                      type="string", metavar="FORMAT",
                      help="Earthquake origin time.",
                      )
    parser.add_option("-c", "--channelcode", dest="channelCode",
                      type="str", metavar="CHANNEL",
                      help="Channel code used by synthetics.",
                      default="HX")
    parser.add_option("-t", "--type", dest="dataType",
                      type="string", metavar="TYPE",
                      help="Type of data to use. [acc,vel,disp]",
                      default="vel")
    parser.add_option("-d", "--datadir", dest="dataDir",
                      type="string", metavar="DIR",
                      help="Directory containing waveform files.",
                      default="OUTPUT_FILES")
    (options, args) = parser.parse_args()
    if len(args) != 0:
        parser.error("Incorrent number of arguments.")

    app = ConvertApp()
    app.filenameIn = options.filenameIn
    app.filenameOut = options.filenameOut
    app.originTime = options.originTime
    app.channelCode = options.channelCode
    app.dataType = options.dataType
    app.dataDir = options.dataDir
    app.run()
  
  
# End of file 
