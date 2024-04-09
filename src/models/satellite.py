import numpy as np
from scipy.integrate import solve_ivp
from .planet import Planet
from shapely.geometry import Polygon
from shapely.ops import cascaded_union


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
        #ja bym dodał altitude, longtitude i to trzecie (wusokość nad planetą, szerokość i długość geograficzna
    
    def update_position(self, new_position, new_velocity):
        self.position = new_position
        self.velocity = new_velocity


    def observe_photo(wysokosc_orbity, szerokosc_geograficzna, dlugosc_geograficzna, szerokosc_pas):
        
        """ 
        Funkcja do obliczania obszaru widzenia satelity Pléiades na Ziemi.
        
        Parametry:
        wysokosc_orbity (float): Wysokość orbity satelity nad Ziemią w kilometrach.
        szerokosc_geograficzna (float): Szerokość geograficzna punktu nad którym znajduje się satelita.
        dlugosc_geograficzna (float): Długość geograficzna punktu nad którym znajduje się satelita.
        szerokosc_pas (float): Szerokość pasma (swath width) przy nadirze w kilometrach.
        
        Zwraca:
        tuple: Kształt i rozmiar obszaru widzenia satelity na Ziemi.
        """
        # Promień Ziemi w kilometrach
        promien_ziemi = 6371
        
        # Obliczenie kąta widzenia na podstawie szerokości pasma i wysokości orbity
        kat_widzenia = 2 * np.arctan((szerokosc_pas / 2) / (wysokosc_orbity + promien_ziemi))
        
        # Obliczenie zasięgu widzenia na powierzchni Ziemi
        zasieg_widzenia = wysokosc_orbity * np.tan(kat_widzenia / 2)
        
        # Obliczenie współrzędnych granicznych obszaru widzenia
        gorna_granica = szerokosc_geograficzna + zasieg_widzenia / 111  # 1 stopień szerokości geograficznej to około 111 km
        dolna_granica = szerokosc_geograficzna - zasieg_widzenia / 111
        prawa_granica = dlugosc_geograficzna + zasieg_widzenia / (111 * np.cos(np.radians(szerokosc_geograficzna)))
        lewa_granica = dlugosc_geograficzna - zasieg_widzenia / (111 * np.cos(np.radians(szerokosc_geograficzna)))
        
        return (gorna_granica, dolna_granica, prawa_granica, lewa_granica)

    def observe_shapely(orbit_altitude, geographic_latitude, geographic_longitude, swath_width):
        # Earth's radius in kilometers
        earth_radius = 6371

        # Calculate the field of view based on the swath width and orbit altitude
        field_of_view = 2 * np.arctan((swath_width / 2) / (orbit_altitude + earth_radius))

        # Calculate the range of view on the Earth's surface
        range_of_view = orbit_altitude * np.tan(field_of_view / 2)

        # Calculate the boundary coordinates of the viewing area
        upper_boundary = geographic_latitude + range_of_view / 111  # 1 degree of geographic latitude is about 111 km
        lower_boundary = geographic_latitude - range_of_view / 111
        right_boundary = geographic_longitude + range_of_view / (111 * np.cos(np.radians(geographic_latitude)))
        left_boundary = geographic_longitude - range_of_view / (111 * np.cos(np.radians(geographic_latitude)))

        # Check if the rectangle crosses the -180/180 degree longitude boundary
        if right_boundary > 180 or left_boundary < -180:
            # Split the rectangle into two parts
            right_rectangle = Polygon(
                [(-180, lower_boundary), (right_boundary % 180, lower_boundary), (right_boundary % 180, upper_boundary),
                 (-180, upper_boundary)])
            left_rectangle = Polygon(
                [(left_boundary % -180, lower_boundary), (180, lower_boundary), (180, upper_boundary),
                 (left_boundary % -180, upper_boundary)])
            return [right_rectangle, left_rectangle]
        else:
            # Return the original rectangle
            return [Polygon(
                [(left_boundary, lower_boundary), (right_boundary, lower_boundary), (right_boundary, upper_boundary),
                 (left_boundary, upper_boundary)])]

    def observe_circle(self, planet_radius):
    
        # Obliczanie odległości do horyzontu z satelity
        odleglosc_do_horyzontu = math.sqrt(self.altitude * (2 * planet_radius + self.altitude))
    
        # Obliczanie promienia obszaru widzenia satelity na powierzchni Ziemi
        promien_obszaru = odleglosc_do_horyzontu * math.tan(math.radians(kat_widzenia / 2))
    
        # Przeliczenie promienia obszaru na stopnie geograficzne
        radius_of_observe = math.degrees(promien_obszaru / planet_radius)
    
        return {
            'srodek': {
                'szerokosc_geograficzna': self.latitude,
                'dlugosc_geograficzna': self.longtitude
            },
            'promien': radius_of_observe
        }
    
    # Przykładowe użycie funkcji
    szerokosc_geograficzna = 51.5074  # szerokość geograficzna Londynu
    dlugosc_geograficzna = -0.1278  # długość geograficzna Londynu
    wysokosc_satelity = 770  # wysokość satelity w kilometrach
    kat_widzenia = 17.6  # kąt widzenia w stopniach
    
    print(obszar_widzenia_satelity(szerokosc_geograficzna, dlugosc_geograficzna, wysokosc_satelity, kat_widzenia))



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

    def calculate_total_area(self):
        total_area = []
        for i in range(len(self.trajectory[0])):
            x, y, z = self.trajectory[0][i], self.trajectory[1][i], self.trajectory[2][i]
            r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
            latitude = np.arcsin(z / r)
            longitude = np.arctan2(y, x)
            area = self.observe_shapely(self.satellite.semi_major_axis - 6371, latitude, longitude,
                                      self.satellite.szerokosc_pas)
            total_area.extend(area)
        return total_area

# Further implementation would follow.
