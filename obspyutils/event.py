# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# ======================================================================


#-----------------------------------------------------------------------
def event_name(event):
    return "".join(event.resource_id.id.split('/')[-2:]).lower()


# ----------------------------------------------------------------------
def find_origin(event, methodId):
    for origin in event.origins:
        if methodId == origin.method_id.id:
            return origin
    raise ValueError("Could not find origin in event %s for method %s." % \
                     (event, methodId))
    return

# ----------------------------------------------------------------------
def moment_tensor(event):
    for fm in event.focal_mechanisms:
        if not fm.moment_tensor.resource_id is None:
            return fm
    raise ValueError("Could not find moment tensor in event %s." % event)
    return

# ----------------------------------------------------------------------
def first_motion(event):
    for fm in event.focal_mechanisms:
        if 'nodal_plane_1' in fm.nodal_planes:
            return fm
    raise ValueError("Could not find focal mechanism with nodal planes in event %s." % event)
    return



# End of file
