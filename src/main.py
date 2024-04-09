from optimization.genetic_algorithm import GeneticAlgorithm
from models.planet import earth, moon

gen_alg = GeneticAlgorithm(planet=earth,
                           constellation_size=50,
                           popSize=20,
                           numParents=8,
                           pm=0.05)

gen_alg.geneticAlgorithm()