class GroundStation:
    def __init__(self, name, latitude, longitude):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude


class Planet:
    def __init__(self, name, radius, mu, rotation_period, j2, angular_velocity, hs_lower_bound, hs_upper_bound):
        self.name = name
        self.radius = radius  # in meters
        self.mu = mu  # Standard gravitational parameter (mu) in m^3/s^2
        self.rotation_period = rotation_period  # Rotation period (sidereal day) in seconds
        self.j2 = j2  # Flattening factor
        self.angular_velocity = angular_velocity  # Angular velocity in rad/s
        self.ground_stations = []
        self.hs_lower_bound = hs_lower_bound  # hs - heliosynchronous
        self.hs_upper_bound = hs_upper_bound

    def add_ground_station(self, name, latitude, longitude):
        self.ground_stations.append(GroundStation(name, latitude, longitude))

    def is_satellite_visible(self, satellite):
        for station in self.ground_stations:
            # Assuming the satellite class has an attribute 'altitude'
            # and a method 'is_visible_from_station' correctly implemented.
            if satellite.is_visible_from_station(station.latitude, station.longitude, satellite.altitude, self.radius):
                return True
        return False


# Example of creating Earth and Moon with their respective constants
earth = Planet(
    name="Earth",
    radius=6371e3,  # in meters
    mu=3.986004418e14,
    rotation_period=86164,
    j2=1.08263e-3,
    angular_velocity=7.2921159e-5,
    hs_lower_bound=6371e3+600,
    hs_upper_bound=6371e3+800,
)

earth.add_ground_station("North", 90, 0)
earth.add_ground_station("South", 90, 0)
earth.add_ground_station("Equator 1", 0, 0)
earth.add_ground_station("Equator 2", 0, 180)

moon = Planet(
    name="Moon",
    radius=17371e2,  # in meters
    mu=4.9048695e12,
    rotation_period=2360591.6,
    j2=202.7e-6,
    angular_velocity=2.6617e-6,
    hs_lower_bound=None,
    hs_upper_bound=None,
)

# Now, both 'earth' and 'moon' objects carry their unique physical constants
# and can be used for simulations specific to each celestial body.
