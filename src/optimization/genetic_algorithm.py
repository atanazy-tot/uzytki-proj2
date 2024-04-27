import math
import random
import numpy as np
from src.models.satellite import Satellite
from src.simulation.coverage import CoverageSimulator


class GeneticAlgorithm:
    def __init__(self,
                 planet,
                 constellation_size: int,
                 popSize: int,
                 numParents: int,
                 pm: float,
                 time_span = 900,
                 swath_width = 800,
                 k: int = 2,
                 numIterations: int = 100):

        assert numParents % 2 == 0

        self.planet = planet
        self.ground_stations = self.planet.ground_stations
        self.constellation_size = constellation_size
        self.popSize = popSize
        self.numParents = numParents
        self.pm = pm  # prob. of mutation
        self.time_span = time_span
        self.k = k  # for parent selection
        self.numIterations = numIterations
        self.minimize = False  # min/max evaluation function
        self.swath_width = swath_width
        self.best_constellation = None

    def random_satellite(self):
        semi_major_axis = random.uniform(self.planet.hs_lower_bound, self.planet.hs_upper_bound)
        eccentricity = random.uniform(0, 0.01)
        inclination = math.radians(98)  # typical for sun-synchronous orbits
        raan = random.uniform(0, 2 * math.pi)
        arg_of_perigee = random.uniform(0, 2 * math.pi)
        true_anomaly = random.uniform(0, 2 * math.pi)
        translation_factor = np.random.randint(0, 5000)
        satellite = Satellite(semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly,
                              self.planet, self.swath_width, translation_factor, self.time_span)
        satellite.update_position(random.choice(range(len(satellite.ground_track))))

        return satellite

    def evaluation(self, constellation):  # In our case, we want to maximize coverage
        self.minimize = False
        c = CoverageSimulator(constellation)
        result = c.calculate_coverage()
        return result

    def mutation(self, constellation):
        for satellite in constellation:
            p = random.random()
            if p < self.pm:  # pm - prob. of mutation
                satellite.semi_major_axis += random.uniform(-10000, 10000)
                satellite.semi_major_axis = max(min(satellite.semi_major_axis, satellite.planet.hs_upper_bound), satellite.planet.hs_lower_bound)
                satellite.eccentricity += random.gauss(0, 0.01)
                satellite.eccentricity = max(min(satellite.eccentricity, 1), 0)  # Eccentricity must be between 0 and 1

                satellite.raan = (satellite.raan + random.gauss(0, 0.01)) % (2 * math.pi)
                satellite.arg_of_perigee = (satellite.arg_of_perigee + random.gauss(0, 0.01)) % (2 * math.pi)
                satellite.true_anomaly = (satellite.true_anomaly + random.gauss(0, 0.01)) % (2 * math.pi)
                satellite.update_position(satellite.translation_factor + random.choice(list(range(-10, 10, 1))))

        return constellation

    def cross_1_point(self, parent_1, parent_2):  # one point crossover
        crossover_point = random.randint(0, len(parent_1))  # parents have the same length (no. satellites)

        child_1 = parent_1[:crossover_point] + parent_2[crossover_point:]
        child_2 = parent_2[:crossover_point] + parent_1[crossover_point:]

        return child_1, child_2

    def cross_bit_flip(self, parent_1, parent_2):  # bit flip crossover
        child_1 = parent_1.copy()
        child_2 = parent_2.copy()

        for satellite in range(self.constellation_size):
            p = random.random()
            if p < 0.5:
                child_1[satellite], child_2[satellite] = child_2[satellite], child_1[satellite]

        return child_1, child_2

    def select_parents(self, fitness):  # tournament selection
        """
        Draws k individuals and returns the best one
        fitness -- fitness matrix for population
        """
        candidates = random.sample(range(self.popSize), self.k)
        best_candidate = min(candidates, key=lambda x: fitness[x]) if self.minimize else max(candidates,
                                                                                             key=lambda x: fitness[x])
        return best_candidate

    def survival(self, population, fitness, offspring, fitness_offspring):
        """
        finds survivors
        """
        population.extend(offspring)
        fitness.extend(fitness_offspring)
        if self.minimize:
            order = sorted(range(len(fitness)), key=lambda k: fitness[k])
        else:
            order = sorted(range(len(fitness)), key=lambda k: -fitness[k])

        # Select the top individuals based on the fitness scores
        survivors = [population[i] for i in order[:self.popSize]]
        fitness_survivors = [fitness[i] for i in order[:self.popSize]]

        return survivors, fitness_survivors

    def geneticAlgorithm(self,
                         cross_1_point: bool = False,
                         stop_crit: bool = True,
                         stop_value: float = 1.0):
        '''
        Genetic algorithm
        stop_crit -- True if the algorithm should be stopped when a certain value is reached
        stop_value -- the value for which the algorithm stops, when stop_crit = True
        '''

        # Initial population
        print("Initializing population...")
        population = [[self.random_satellite() for _ in range(self.constellation_size)] for _ in range(self.popSize)]
        print("Calculating fitness...")
        fitness = [self.evaluation(constellation) for constellation in population]
        print("Initial:", "min = ", min(fitness), ", max = ", max(fitness), ", average = ", sum(fitness)/len(fitness))

        # Main loop. Stop criterion: number of iterations
        for num_it in range(self.numIterations):
            if stop_crit and (
                    (self.minimize and min(fitness) == stop_value) or (not self.minimize and max(fitness) == stop_value)):
                break

            offspring = []

            for kk in range(self.numParents // 2):
                # Parent selection
                parent_1 = population[self.select_parents(fitness)]
                parent_2 = population[self.select_parents(fitness)]

                # Crossover
                if cross_1_point:
                    child_1, child_2 = self.cross_1_point(parent_1, parent_2)
                else:
                    child_1, child_2 = self.cross_bit_flip(parent_1, parent_2)

                # Mutation
                child_1 = self.mutation(child_1)
                child_2 = self.mutation(child_2)

                # Save the children
                offspring.extend([child_1, child_2])

            # Evaluation of children
            fitness_offspring = [self.evaluation(constellation) for constellation in offspring]

            # Survivor selection
            population, fitness = self.survival(population, fitness, offspring, fitness_offspring)

            # Print info
            print("Iteration ", num_it, ": min = ", min(fitness), ", max = ", max(fitness), ", average = ",
                  sum(fitness)/len(fitness))

            arg = fitness.index(min(fitness)) if self.minimize else fitness.index(max(fitness))
            self.best_constellation = population[arg]

            # for printing orbit parameters
            # Open the file in write mode to erase its contents
            with open("sat_params.txt", "w") as f:
                pass  # Do nothing

            # Now open the file in append mode to add new data
            with open("sat_params.txt", "a") as f:
                for s in self.best_constellation:
                    f.write(str(s) + '\n')  # Convert to string and append

        arg = fitness.index(min(fitness)) if self.minimize else fitness.index(max(fitness))
        return population[arg]

