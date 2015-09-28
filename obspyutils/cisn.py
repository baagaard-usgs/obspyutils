# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================


#-----------------------------------------------------------------------
def remove_structures(inventory):
    """
    Remove structures from inventory of CISN stations using daily updated
    adhoc channel list.

    :param: inventory Inventory of stations.
    """

    import urllib2
    adhocURL = "http://agent86.gps.caltech.edu/utils/adhoc.lis"
    fin = urllib2.urlopen(adhocURL)
    lines = fin.readlines()
    fin.close()
    
    structures = []
    for line in lines:
        if line[0] == "#" or len(line.strip()) == 0:
            continue
        fields = line.split()
        if len(fields) < 9:
            raise ValueError("Could not split line '%s'." % line.strip())
        station = fields[0]
        network = fields[1]
        channel = fields[2]
        location = fields[3]
        cosmos = int(fields[4])
        lat = float(fields[5])
        lon = float(fields[6])
        elev = float(fields[7])
        description = fields[8]
        if cosmos >= 5:
            key = "%s.%s" % (network, station)
            if not key in structures:
                structures.append(key)

    for network in inventory.networks:
        remove = []
        for station in network.stations:
            key = "%s.%s" % (network.code, station.code)
            if key in structures:
                remove.append(station)
        for station in remove:
            network.stations.remove(station)


if __name__ == "__main__":
    import obspy.core.inventory
    inventory = obspy.core.inventory.Inventory([], "")
    remove_structures(inventory)

# End of file
