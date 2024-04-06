class GroundStation:
    def __init__(self, name, latitude, longitude):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude


class Planet:
    def __init__(self, name, radius):
        self.name = name
        self.radius = radius
        self.ground_stations = []

    def add_ground_station(self, name, latitude, longitude):
        self.ground_stations.append(GroundStation(name, latitude, longitude))

    def is_satellite_visible(self, satellite):
        for station in self.ground_stations:
            if satellite.is_visible_from_station(station.latitude, station.longitude, self.radius):
                return True
        return False
