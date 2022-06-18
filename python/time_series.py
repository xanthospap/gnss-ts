from car2ell import cartesian2ellipsoidal

class TimeSeries:

    def __init__(self, site, entries=[]):
        self._site = site
        ## list of dictionaries ...
        self._data = entries

    def sort(self, onkey='t'):
        ## onkey can be any of the keys included in the self._data list 
        ## elements
        return TimeSeries(self._site, sorted(self._data, key=lambda d return d[onkey]))

    def average_crd(self):
        ax = ay = az = 0e0
        n = 0
        for entry in self._data:
            if not entry['ignore']:
                ax += (entry['x'] - ax) / (n+1)
                ay += (entry['y'] - ay) / (n+1)
                az += (entry['z'] - az) / (n+1)
        return ax,ay,az

    def topocentric(self):
        x, y, z = self.average_crd()
        lon, lat, hgt = cartesian2ellipsoidal(x,y,z)

        nd = [ {cartesian2ellipsoidal(entry['x'], entry['y'], entry['z'])} for entry in self._data ]
