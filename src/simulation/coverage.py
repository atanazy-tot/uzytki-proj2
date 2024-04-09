class CoverageSimulator:
    def __init__(self, satellites, ground_stations, planet):
        self.satellites = satellites
        self.ground_stations = ground_stations
        self.planet = planet

    def calculate_coverage(self, start_time, end_time, time_step):
        tot_area = []
        for satellite in self.satelites:
            p = OrbitPropagator(satellite, start_time, end_time, time_step)
            p.propagate_orbit()
            area = p.calculate_total_area()
            tot_area.extend(area)

        total_area_polygon = cascaded_union(tot_area)

        # Calculate the total area covered by the photos
        total_area_covered = total_area_polygon.area

        # Calculate the proportion covered
        whole_map = Polygon([(-180, -90), (180, -90), (180, 90), (-180, 90)])
        map_area = whole_map.area
        proportion_covered = total_area_covered / earth_surface_area

        return proportion_covered