from shapely.ops import cascaded_union
from shapely.geometry import Polygon
from shapely import unary_union

class CoverageSimulator:
    def __init__(self, constellation):
        self.constellation = constellation

    def calculate_coverage(self):
        tot_area = []
        for satellite in self.constellation:
            for i in range(len(satellite.sat_trajectory)):
                latitude, longitude, altitude = satellite.sat_trajectory[i]
                area = satellite.observe_shapely(altitude, latitude, longitude)
                for polygon in area:
                    tot_area.append(polygon)
        total_area_polygon = unary_union(tot_area)

        # Calculate the total area covered by the photos
        total_area_covered = total_area_polygon.area

        # Calculate the proportion covered
        whole_map = Polygon([(-180, -90), (180, -90), (180, 90), (-180, 90)])
        map_area = whole_map.area
        proportion_covered = total_area_covered / map_area

        return proportion_covered
