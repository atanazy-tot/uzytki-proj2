import numpy as np


def spherical_to_cartesian(r, theta, phi):
    # zakładamy wartości kątów w radianach
    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)
    return np.array([x, y, z])


def cartesian_to_spherical(x, y, z):
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y, x)
    phi = np.arcsin(z / r)
    return np.array([r, theta, phi])


def geographic_to_eci(latitude_deg, longitude_deg, altitude):
    pass


def eci_to_geographic():
    pass