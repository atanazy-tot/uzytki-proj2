from .planet import Planet


class Satellite:
    def __init__(self,
                 name: str,
                 mass: float,
                 planet: Planet,
                 latitude: float,
                 longitude: float,
                 altitude: float,
                 velocity: float,
                 angle: float):
        self.name = name
        self.mass = mass
        self.planet = planet
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude
        self.velocity = velocity  # prędkość początkowa? nie jest wyznaczana przez heliosynchroniczność?
        self.angle = angle  # czy to nie jest jedyne, co nam wystarczy przy heliosynchroniczności?

    def update_position(self, time_delta):
        orbital_circumference = 2 * 3.14159 * (self.planet.radius + self.altitude)
        orbit_period = orbital_circumference / self.velocity
        degrees_per_second = 360 / orbit_period
        degrees_moved = degrees_per_second * time_delta

        self.longitude += degrees_moved
        self.longitude %= 360
