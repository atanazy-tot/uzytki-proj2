import numpy as np


def calculate_orbit_parameters(latitude_deg, longitude_deg, altitude_m, local_time_h):
    # Constants
    mu = 398600.4418  # Earth's gravitational parameter, km^3/s^2
    R_earth = 6371  # Earth's radius, km
    J2 = 1.08263e-3  # Earth's J2 value
    inclination = np.radians(98.6)  # Typical SSO inclination in radians

    # Convert altitude from meters to kilometers
    altitude_km = altitude_m / 1000

    # Calculate semi-major axis
    semi_major_axis = R_earth + altitude_km

    # Orbital period (in seconds)
    T = 2 * np.pi * np.sqrt(semi_major_axis ** 3 / mu)

    # Approximate local solar time to RAAN conversion
    # This is a simplified model and doesn't account for the exact dynamics
    # Assuming the satellite crosses the equator at the local noon (maximum simplification)
    GMT_offset = longitude_deg / 15  # Convert longitude to GMT offset
    RAAN_offset_hours = local_time_h - 12  # Difference from local noon
    daily_precession_deg = 360 / 365.25  # Daily precession needed for SSO
    RAAN = (GMT_offset + RAAN_offset_hours) * 15 - daily_precession_deg  # Simplified RAAN calculation

    # Ensure RAAN is within 0 to 360 degrees
    RAAN = RAAN % 360

    return {
        "semi_major_axis_km": semi_major_axis,
        "inclination_deg": np.degrees(inclination),
        "RAAN_deg": RAAN,
        "period_s": T
    }


# Example usage: Calculate parameters for a satellite to pass over Warsaw at 10:00 local time
latitude_deg = 52.2297
longitude_deg = 21.0122
altitude_m = 700000  # 700 km altitude
local_time_h = 10.0  # 10:00 AM

orbit_params = calculate_orbit_parameters(latitude_deg, longitude_deg, altitude_m, local_time_h)
print(orbit_params)
