from optimization.genetic_algorithm import GeneticAlgorithm
from models.planet import earth, moon

earth.add_ground_station("North", 90, 0)
earth.add_ground_station("South", 90, 0)
earth.add_ground_station("Equator 1", 0, 0)
earth.add_ground_station("Equator 2", 0, 180)

gen_alg = GeneticAlgorithm(planet=earth,
                           constellation_size=8,
                           popSize=10,
                           numParents=8,
                           pm=0.05)

gen_alg.geneticAlgorithm()

