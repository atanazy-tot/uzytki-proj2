class GroundStation:
    def __init__(self, name, latitude, longitude):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude


class Planet:
    def __init__(self, name, radius, mu, rotation_period, j2, angular_velocity):
        self.name = name
        self.radius = radius  # in meters
        self.mu = mu  # Standard gravitational parameter (mu) in m^3/s^2
        self.rotation_period = rotation_period  # Rotation period (sidereal day) in seconds
        self.j2 = j2  # Flattening factor
        self.angular_velocity = angular_velocity  # Angular velocity in rad/s
        self.ground_stations = []

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
    radius=6371e3,
    mu=3.986004418e14,
    rotation_period=86164,
    j2=1.08263e-3,
    angular_velocity=7.2921159e-5
)

moon = Planet(
    name="Moon",
    radius=17371e2,
    mu=4.9048695e12,
    rotation_period=2360591.6,
    j2=202.7e-6,
    angular_velocity=2.6617e-6
)

# Now, both 'earth' and 'moon' objects carry their unique physical constants
# and can be used for simulations specific to each celestial body.
