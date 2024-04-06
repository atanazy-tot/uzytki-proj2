class Planet:
    def __init__(self, gravitational_constant, radius, j2):
        """
        Initializes a planet object with the necessary gravitational parameters.
        :param gravitational_constant: Gravitational constant (G * M) for the planet in m^3/s^2
        :param radius: Mean radius of the planet in meters
        :param j2: Second zonal harmonic coefficient for the planet's gravitational field
        """
        self.gravitational_constant = gravitational_constant
        self.radius = radius
        self.j2 = j2
