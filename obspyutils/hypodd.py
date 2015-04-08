# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================

import obspy

#-----------------------------------------------------------------------
def tocatalog(filename="tomoDD.reloc"):
    """
    Collect HypoDD/tomoDD output into a obspy catalog.
    """
    fin = open(filename, "r")
    lines = fin.readlines()
    fin.close()

    events = []
    for line in lines:
        fields = line.split()
        if len(fields) != 24:
            import pdb
            pdb.set_trace()
            raise IOError("Unrecognized format for HypoDD/tomoDD file.\nLine: '%s'" % line)
        
        evid = fields[0].strip()
        latitude = float(fields[1])
        longitude = float(fields[2])
        depth = float(fields[3])
        
        x = float(fields[4])
        y = float(fields[5])
        z = float(fields[6])
        ex = float(fields[7])
        ey = float(fields[8])
        ez = float(fields[9])

        # Date/time is 10-15
        year = int(fields[10])
        month = int(fields[11])
        day = int(fields[12])
        hour = int(fields[13])
        minute = int(fields[14])
        second = int(fields[15])/100
        fracsecond = int(fields[15]) % 100
        datetime = obspy.UTCDateTime("%d-%d-%dT%d:%d:%d.%d" % (year, month, day, hour, minute, second, fracsecond))

        # mag is 16
        magnitude = float(fields[16])
        
        # Unknown fields 17,18,19,20,21,22
        nccp = int(fields[17])
        nccs = int(fields[18])
        nctp = int(fields[19])
        ncts = int(fields[20])
        rcc = float(fields[21])
        rct = float(fields[22])
        cid = int(fields[23])

        event = obspy.core.event.Event()
        
        origin = obspy.core.event.Origin()
        origin.time = datetime
        origin.latitude = latitude
        origin.longitude = longitude
        origin.depth = depth * 1000
        origin.depth_type = "from location"
        origin.method_id = obspy.core.event.ResourceIdentifier(id="smi:hardebeck.usgs.gov/origin/TOMODD")
        origin.type = "hypocenter"
        event.origins.append(origin)

        events.append(event)
                
    catalog = obspy.core.event.Catalog(events=events)
    return catalog




# End of file
