import numpy as np
from scipy.integrate import solve_ivp
from .planet import Planet


class Orbit:
    def __init__(self, semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly):
        self.semi_major_axis = semi_major_axis  # in kilometers
        self.eccentricity = eccentricity
        self.inclination = inclination  # in radians
        self.raan = raan  # Right Ascension of the Ascending Node, in radians
        self.arg_of_perigee = arg_of_perigee  # in radians
        self.true_anomaly = true_anomaly  # in radians

    def to_state_vector(self):
        # Convert the orbital elements to a state vector (position and velocity)
        def to_state_vector(self):
            # Constants
            mu = 398600.4418  # Earth’s gravitational parameter, km^3/s^2

            # Position in perifocal coordinates
            p = self.semi_major_axis * (1 - self.eccentricity ** 2)
            r = p / (1 + self.eccentricity * np.cos(self.true_anomaly))
            r_pqw = np.array([r * np.cos(self.true_anomaly),
                              r * np.sin(self.true_anomaly),
                              0])

            # Velocity in perifocal coordinates
            v_pqw = np.array([-np.sqrt(mu / p) * np.sin(self.true_anomaly),
                              np.sqrt(mu / p) * (self.eccentricity + np.cos(self.true_anomaly)),
                              0])

            # Rotation matrices
            R_3_O = np.array([[np.cos(-self.raan), -np.sin(-self.raan), 0],
                              [np.sin(-self.raan), np.cos(-self.raan), 0],
                              [0, 0, 1]])

            R_1_i = np.array([[1, 0, 0],
                              [0, np.cos(-self.inclination), -np.sin(-self.inclination)],
                              [0, np.sin(-self.inclination), np.cos(-self.inclination)]])

            R_3_w = np.array([[np.cos(-self.arg_of_perigee), -np.sin(-self.arg_of_perigee), 0],
                              [np.sin(-self.arg_of_perigee), np.cos(-self.arg_of_perigee), 0],
                              [0, 0, 1]])

            # Perform the rotations to get into the ECI frame
            r_eci = R_3_O @ R_1_i @ R_3_w @ r_pqw
            v_eci = R_3_O @ R_1_i @ R_3_w @ v_pqw

            return r_eci, v_eci
        pass


class Satellite(Orbit):
    def __init__(self, semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly,
                 resolution, optical_bands, camera_fov, planet: Planet):
        super().__init__(semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly)
        self.resolution = resolution
        self.optical_bands = optical_bands
        self.camera_fov = camera_fov
        self.planet = planet
        self.position = None
        self.velocity = None

    def update_position(self, new_position, new_velocity):
        self.position = new_position
        self.velocity = new_velocity

    def observe(self):
        # Simulate the satellite taking an observation with its camera
        pass

    def is_visible_from_station(self, satellite_lat, satellite_long, satellite_alt, planet_radius):
        # Method to determine visibility from a ground station
        # Obliczanie odległości między odbiornikiem a satelitą
        d_lat = math.radians(satellite_lat - receiver_lat)
        d_lon = math.radians(satellite_lon - receiver_lon)
        a = math.sin(d_lat/2) * math.sin(d_lat/2) + math.cos(math.radians(receiver_lat)) * math.cos(math.radians(satellite_lat)) * math.sin(d_lon/2) * math.sin(d_lon/2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        d = planet_radius * c
    
        # Obliczanie kąta między odbiornikiem a satelitą
        angle_to_satellite = math.atan2(satellite_alt - receiver_alt, d)
    
        # Obliczanie kąta między odbiornikiem a horyzontem
        angle_to_horizon = math.acos(planet_radius / (planet_radius + receiver_alt))

        # Sprawdzanie, czy satelita jest widoczny
        return angle_to_satellite < angle_to_horizon

class OrbitPropagator:
    def __init__(self, satellite, start_time, end_time, time_step):
        self.satellite = satellite
        self.start_time = start_time
        self.end_time = end_time
        self.time_step = time_step
        self.trajectory = []

    def propagate_orbit(self):
        # Define the differential equations for the two-body problem
        def equations_of_motion(t, y):
            # Constants
            mu = 398600.4418  # Earth’s gravitational parameter, km^3/s^2
            J2 = 1.08263e-3
            R_earth = 6378.137  # Earth’s equatorial radius, km

            # Unpack the position and velocity vectors
            r_vec = y[:3]
            v_vec = y[3:]
            r = np.linalg.norm(r_vec)

            # Two-body acceleration
            acc_gravity = -mu * r_vec / r**3

            # J2 perturbation
            z2 = r_vec[2]**2
            r2 = r**2
            tx = r_vec[0] / r * (5 * z2 / r2 - 1)
            ty = r_vec[1] / r * (5 * z2 / r2 - 1)
            tz = r_vec[2] / r * (5 * z2 / r2 - 3)
            acc_j2 = 1.5 * J2 * mu * R_earth**2 / r**4 * np.array([tx, ty, tz])

            # Total acceleration
            acc_total = acc_gravity + acc_j2

            return np.hstack((v_vec, acc_total))

        # Set the initial conditions for the integrator
        y0 = np.hstack(self.satellite.to_state_vector())
        t_span = (self.start_time, self.end_time)
        t_eval = np.arange(self.start_time, self.end_time, self.time_step)

        # Solve the system of differential equations
        solution = solve_ivp(equations_of_motion, t_span, y0, t_eval=t_eval, rtol=1e-6)

        # Store the results in the trajectory attribute
        self.trajectory = solution.y
        return solution.y

    def calculate_ground_track(self):
        # Earth's equatorial radius in kilometers
        R_earth = 6378.137
        # Earth's rotation rate (radians per second)
        omega_earth = 7.2921159e-5

        ground_track = []

        for i in range(len(self.trajectory[0])):
            x, y, z = self.trajectory[0][i], self.trajectory[1][i], self.trajectory[2][i]

            # Calculate the current longitude, taking into account Earth's rotation
            time = self.start_time + i * self.time_step
            theta = (omega_earth * time) % (2 * np.pi)
            longitude = np.arctan2(y, x) - theta
            # Ensure longitude is between -pi and pi
            longitude = (longitude + np.pi) % (2 * np.pi) - np.pi

            # Calculate latitude from the position
            r = np.sqrt(x**2 + y**2 + z**2)
            latitude = np.arcsin(z / r)

            # Convert latitude and longitude to degrees
            latitude_deg = np.degrees(latitude)
            longitude_deg = np.degrees(longitude)

            ground_track.append((latitude_deg, longitude_deg))

        return ground_track

# Further implementation would follow.
