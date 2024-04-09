import numpy as np
from scipy.integrate import solve_ivp
from .planet import Planet


class Orbit:
    def __init__(self, semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly,
                 planet: Planet):
        self.semi_major_axis = semi_major_axis  # in kilometers
        self.eccentricity = eccentricity
        self.inclination = inclination  # in radians
        self.raan = raan  # Right Ascension of the Ascending Node, in radians
        self.arg_of_perigee = arg_of_perigee  # in radians
        self.true_anomaly = true_anomaly  # in radians
        self.planet = planet

    def to_state_vector(self):
        # Constants
        mu = self.planet.mu  # Earth’s gravitational parameter, m^3/s^2

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



class Satellite(Orbit):
    def __init__(self, semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly,
                 resolution, optical_bands, camera_fov, swath_width, planet: Planet):
        super().__init__(semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly,
                         planet)
        self.resolution = resolution
        self.optical_bands = optical_bands
        self.swath_width = swath_width
        self.camera_fov = camera_fov
        self.position = None
        self.velocity = None
        self.latitude = None
        self.longitude = None
        self.altitude = None
    
    def update_position(self, latitude, longitude, altitude):
        self.latitude = latitude
        self.longitude = longitude
        self.altitude = altitude

    def observe_area(self, orbit_height, geographic_latitude, geographic_longitude, swath_width):
        """
        Function to calculate the visibility area of a Pleiades satellite on Earth.

        Parameters:
        orbit_height (float): Orbit height of the satellite above Earth in kilometers.
        geographic_latitude (float): Geographic latitude of the point directly below the satellite.
        geographic_longitude (float): Geographic longitude of the point directly below the satellite.
        swath_width (float): Width of the swath at nadir in kilometers.

        Returns:
        tuple: Shape and size of the satellite's visibility area on Earth.
        """
        # Earth's radius in kilometers
        earth_radius = 6371

        # Calculate the viewing angle based on the swath width and orbit height
        viewing_angle = 2 * np.arctan((swath_width / 2) / (orbit_height + earth_radius))

        # Calculate the visibility range on the Earth's surface
        visibility_range = orbit_height * np.tan(viewing_angle / 2)

        # Calculate the coordinates of the visibility area boundaries
        upper_boundary = geographic_latitude + visibility_range / 111  # 1 degree of latitude is approximately 111 km
        lower_boundary = geographic_latitude - visibility_range / 111
        right_boundary = geographic_longitude + visibility_range / (111 * np.cos(np.radians(geographic_latitude)))
        left_boundary = geographic_longitude - visibility_range / (111 * np.cos(np.radians(geographic_latitude)))

        return (upper_boundary, lower_boundary, right_boundary, left_boundary)

    def observe_circle(self, planet_radius, viewing_angle_deg):
        """
        Calculate the observation area as a circle on the planet surface.

        Parameters:
        planet_radius (float): Radius of the planet.
        viewing_angle_deg (float): Viewing angle in degrees.

        Returns:
        dict: Center and radius of the observation area.
        """
        # Calculating distance to the horizon from the satellite
        distance_to_horizon = np.sqrt(self.altitude * (2 * planet_radius + self.altitude))

        # Calculating the radius of the visible area on the planet surface
        visible_area_radius = distance_to_horizon * np.tan(np.radians(viewing_angle_deg / 2))

        # Convert the area radius to geographic degrees
        radius_in_degrees = np.degrees(visible_area_radius / planet_radius)

        return {
            'center': {
                'latitude': self.latitude,
                'longitude': self.longitude
            },
            'radius': radius_in_degrees
        }


    def is_visible_from_station(self, satellite_lat, satellite_long, satellite_alt, receiver_lat, receiver_lon,
                                receiver_alt, planet_radius):
        """
        Determine if the satellite is visible from a ground station.

        Parameters:
        satellite_lat (float): Satellite's latitude.
        satellite_long (float): Satellite's longitude.
        satellite_alt (float): Satellite's altitude.
        receiver_lat (float): Receiver's (ground station) latitude.
        receiver_lon (float): Receiver's (ground station) longitude.
        receiver_alt (float): Receiver's (ground station) altitude.
        planet_radius (float): Radius of the planet.

        Returns:
        bool: True if the satellite is visible from the ground station, False otherwise.
        """
        # Calculating the distance between the receiver and the satellite
        d_lat = np.radians(satellite_lat - receiver_lat)
        d_lon = np.radians(satellite_long - receiver_lon)
        a = np.sin(d_lat / 2) ** 2 + np.cos(np.radians(receiver_lat)) * np.cos(np.radians(satellite_lat)) * np.sin(
            d_lon / 2) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        distance = planet_radius * c

        # Calculating the angle to the satellite
        angle_to_satellite = np.arctan2(satellite_alt - receiver_alt, distance)

        # Calculating the angle to the horizon
        angle_to_horizon = np.arccos(planet_radius / (planet_radius + receiver_alt))

        # Checking if the satellite is visible
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
            mu = self.satellite.planet.mu  # Earth’s gravitational parameter, m^3/s^2
            J2 = self.satellite.planet.J2
            radius = self.satellite.planet.radius  # Earth’s equatorial radius, m

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
            acc_j2 = 1.5 * J2 * mu * radius**2 / r**4 * np.array([tx, ty, tz])

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
        # Earth's rotation rate (radians per second)
        omega_earth = self.satellite.planet.angular_velocity

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
            r = np.sqrt(x**2 + y**2 + z**2)  # that's just altitude bro
            latitude = np.arcsin(z / r)

            # Convert latitude and longitude to degrees
            latitude_deg = np.degrees(latitude)
            longitude_deg = np.degrees(longitude)

            ground_track.append((latitude_deg, longitude_deg, r))

        return ground_track

# Further implementation would follow.
