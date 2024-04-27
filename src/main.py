from optimization.genetic_algorithm import GeneticAlgorithm
from models.planet import earth

gen_alg = GeneticAlgorithm(planet=earth,
                           constellation_size=30,
                           popSize=20,
                           numParents=8,
                           pm=0.3,
                           time_span=10800,
                           swath_width=800,
                           numIterations=25,
)

gen_alg.geneticAlgorithm(cross_1_point=False, file_name = "sat_params_30.txt")  # use bit flip

# Print orbital parameters for the best constellation
print("final constellation:\n")

for i, satellite in enumerate(gen_alg.best_constellation):
    print(i, satellite)
