import math

def topocentric_matrix(lon, lat):
    cf = math.cos(lat)
    sf = math.sin(lat)
    cl = math.cos(lon)
    sl = math.sin(lon)

    ## unit vector along east
    d00 = -sl
    d01 = cl
    d02 = 0e0
    ## unit vector along north
    d10 = -sf * cl
    d11 = -sf * sl
    d12 = cf
    ## unit vector along up
    d20 = cf*cl
    d21 = cf*sl
    d22 = sf

    return  [[d00,d01,d02],[d10,d11,d12],[d20,d21,d22]]
