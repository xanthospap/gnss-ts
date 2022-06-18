import math

def cartesian2ellipsoidal(x, y, z, a=6378137.0e0, f=1e0/298.257223563e0):

  ## Functions of ellipsoid parameters.
  aeps2 = a*a*1e-32
  e2    = (2.0e0-f)*f
  e4t   = e2*e2*1.5e0
  ep2   = 1.0e0-e2
  ep    = math.sqrt(ep2)
  aep   = a*ep

  ## Compute distance from polar axis squared.
  p2 = x*x + y*y

  ## Compute longitude lon
  if p2 != .0e0:
    lon = math.atan2(y, x)
  else:
    lon = .0e0

  ## Ensure that Z-coordinate is unsigned.
  absz = abs(z)

  if p2 > aeps2: ## Continue unless at the poles

    ## Compute distance from polar axis.
    p = math.sqrt(p2)
    ## Normalize.
    s0 = absz/a
    pn = p/a
    zp = ep*s0 
    ## Prepare Newton correction factors.
    c0  = ep*pn 
    c02 = c0*c0 
    c03 = c02*c0 
    s02 = s0*s0 
    s03 = s02*s0 
    a02 = c02+s02 
    a0  = math.sqrt(a02) 
    a03 = a02*a0 
    d0  = zp*a03 + e2*s03 
    f0  = pn*a03 - e2*c03 
    ## Prepare Halley correction factor.
    b0 = e4t*s02*c02*pn*(a0-ep) 
    s1 = d0*f0 - b0*s0 
    cp = ep*(f0*f0-b0*c0) 
    ## Evaluate latitude and height.
    lat = math.atan(s1/cp)
    s12 = s1*s1 
    cp2 = cp*cp 
    hgt = (p*cp+absz*s1-a*math.sqrt(ep2*s12+cp2))/math.sqrt(s12+cp2);

  else: ## Special case: pole.

    lat = math.pi / 2e0;
    hgt = absz - aep;

  ## Restore sign of latitude.
  if z < 0.e0:
    lat = -lat

  return lon, lat, hgt
