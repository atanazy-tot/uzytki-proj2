import numpy as np


def calculate_inclination(n: np.array, z=np.array([0, 0, 1])):
    """
    Oblicza inklinację orbity na podstawie wektora normalnego do płaszczyzny orbity (n)
    oraz wektora osi Z (z), który domyślnie jest ustawiony jako [0, 0, 1].

    :param n: Wektor normalny do płaszczyzny orbity, np.array([n_x, n_y, n_z]).
    :param z: Wektor osi Z, domyślnie ustawiony na osi obrotu planety, np.array([0, 0, 1]).
    :return: Inklinacja orbity w stopniach.
    """
    norm_n = np.linalg.norm(n)
    norm_z = np.linalg.norm(z)

    # Obliczenie cosinusa kąta inklinacji
    cos_theta = np.dot(n, z) / (norm_n * norm_z)

    # Obliczenie inklinacji jako arcus cosinus cos_theta
    inclination_radians = np.arccos(cos_theta)

    # Konwersja na stopnie
    inclination_degrees = np.degrees(inclination_radians)

    return inclination_degrees



