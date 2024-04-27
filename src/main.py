from optimization.genetic_algorithm import GeneticAlgorithm
from models.planet import earth, moon

earth.add_ground_station("North", 90, 0)
earth.add_ground_station("South", 90, 0)
earth.add_ground_station("Equator 1", 0, 0)
earth.add_ground_station("Equator 2", 0, 180)

gen_alg = GeneticAlgorithm(planet=earth,
                           constellation_size=4,
                           popSize=4,
                           numParents=8,
                           pm=0.3,
                           time_span=900,
                           swath_width=800,
                           numIterations=3,
)

gen_alg.geneticAlgorithm(cross_1_point=False)  # use bit flip
print("final constellation:\n")

for i, satellite in enumerate(gen_alg.best_constellation):
    print(i, satellite)
