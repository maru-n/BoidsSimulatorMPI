################################################
# submodule to calculate distance in periodic boundary for DBSCAN
################################################


def periodic_boundary(X, Y, field_size=1.0):
    cdef float xd = abs(X[0] - Y[0])
    xd = min(xd, field_size-xd)
    cdef float yd = abs(X[1] - Y[1])
    yd = min(yd, field_size-yd)
    cdef float zd = abs(X[2] - Y[2])
    zd = min(zd, field_size-zd)
    cdef float d = (xd*xd + yd*yd + zd*zd)**(0.5)
    return d


def calc_distance_xv(X, Y, alpha=0.5):
    cdef float xd = abs(X[0] - Y[0])
    xd = min(xd, 1.0-xd)
    cdef float yd = abs(X[1] - Y[1])
    yd = min(yd, 1.0-yd)
    cdef float zd = abs(X[2] - Y[2])
    zd = min(zd, 1.0-zd)
    vxd = (X[3] - Y[3])
    vyd = (X[4] - Y[4])
    vzd = (X[5] - Y[5])
    cdef float d = alpha * (xd*xd + yd*yd + zd*zd)**(0.5) + (1.0-alpha) * (vxd*vxd + vyd*vyd + vzd*vzd)**(0.5)
    return d

def calc_distance_xs(X, Y, alpha=0.5):
    cdef float xd = abs(X[0] - Y[0])
    xd = alpha * min(xd, 1.0-xd)
    cdef float yd = abs(X[1] - Y[1])
    yd = alpha * min(yd, 1.0-yd)
    cdef float zd = abs(X[2] - Y[2])
    zd = alpha * min(zd, 1.0-zd)
    cdef float sd = (1.0-alpha) * abs(X[3] - Y[3])
    #cdef float d = alpha * (xd*xd + yd*yd + zd*zd)**(0.5) + (1.0-alpha) * sd
    cdef float d = (xd*xd + yd*yd + zd*zd + sd*sd)**(0.5)
    return d

def calc_distance_xs_legacy(X, Y, field_size=1.0, alpha=0.5):
    cdef float xd = abs(X[0] - Y[0])
    xd = min(xd, field_size-xd)
    cdef float yd = abs(X[1] - Y[1])
    yd = min(yd, field_size-yd)
    cdef float zd = abs(X[2] - Y[2])
    zd = min(zd, field_size-zd)
    sd = X[3] - Y[3]
    cdef float d = alpha * (xd*xd + yd*yd + zd*zd)**(0.5) + (1.0-alpha) * sd
    return d
