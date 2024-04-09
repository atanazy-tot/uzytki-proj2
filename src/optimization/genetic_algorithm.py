import math
import random
import numpy as np

class GeneticAlgorithm:
    def __init__(self,
                 ground_stations,
                 planet,
                 start_time,
                 end_time,
                 time_step,
                 constellation_size: int,
                 popSize: int,
                 numParents: int,
                 pm: float,
                 k: int = 2,
                 numIterations: int = 100):

        assert numParents % 2 == 0

        self.ground_stations = ground_stations
        self.planet = planet
        self.start_time = start_time
        self.end_time = end_time
        self.time_step = time_step
        self.constellation_size = constellation_size
        self.popSize = popSize
        self.numParents = numParents
        self.pm = pm # prob. of mutation
        self.k = k # for parent selection
        self.numIterations = numIterations
        self.minimize = False # min/max evaluation function

    def random_satellite(self):
        semi_major_axis = random.uniform(self.planet.hs_lower_bound, self.planet.hs_upper_bound)
        eccentricity = random.uniform(0, 0.01)
        inclination = math.radians(98)  # typical for sun-synchronous orbits
        raan = random.uniform(0, 2 * math.pi)
        arg_of_perigee = random.uniform(0, 2 * math.pi)
        true_anomaly = random.uniform(0, 2 * math.pi)

        return Satellite(semi_major_axis, eccentricity, inclination, raan, arg_of_perigee, true_anomaly, self.planet)

    def evaluation(self, constellation): # In our case, we want to maximize coverage
        self.minimize = False
        c = CoverageSimulator(constellation, self.ground_stations, self.planet)
        result = c.calculate_covarage(self.start_time, self.end_time, self.time_step)
        return result

    def mutation(self, constellation):
        for satellite in constellation:
            p = random.random()
            if p < self.pm: # pm - prob. of mutation
                satellite.semi_major_axis += random.uniform(-10, 10)
                satellite.semi_major_axis = max(min(satellite.semi_major_axis, satellite.planet.hs_upper_bound), satellite.planet.hs_lower_bound)
                satellite.eccentricity += random.gauss(0, 0.01)
                satellite.eccentricity = max(min(satellite.eccentricity, 1), 0)  # Eccentricity must be between 0 and 1

                satellite.raan = (satellite.raan + random.gauss(0, 0.01)) % (2 * math.pi)
                satellite.arg_of_perigee = (satellite.arg_of_perigee + random.gauss(0, 0.01)) % (2 * math.pi)
                satellite.true_anomaly = (satellite.true_anomaly + random.gauss(0, 0.01)) % (2 * math.pi)

        return constellation

    def cross(self, parent_1, parent_2): # one point crossover
        crossover_point = random.randint(0, len(parent_1)) # parents have the same length (no. satellites)

        child_1 = parent_1[:crossover_point] + parent_2[crossover_point:]
        child_2 = parent_2[:crossover_point] + parent_1[crossover_point:]

        return child_1, child_2

    def select_parents(self, fitness): # tournament selection
        """
        Draws k individuals and returns the best one
        fitness -- fitness matrix for population
        """
        candidates = random.sample(range(self.popSize), self.k)
        best_candidate = np.min(candidates, key=lambda x: fitness[x]) if self.minimize else max(candidates,
                                                                                             key=lambda x: fitness[x])
        return best_candidate

    def survival(self, population, fitness, offspring, fitness_offspring):
        """
        finds survivors
        """
        current_generation = population.extend(offspring)
        fitness_generation = fitness.extend(fitness_offspring)

        order = sorted(range(len(fitness_generation)), key=lambda k: fitness_generation[k])

        # Select the top individuals based on the fitness scores
        survivors = [current_generation[i] for i in order[:self.popSize]]
        fitness_survivors = [fitness_generation[i] for i in order[:self.popSize]]

        return survivors, fitness_survivors

    def geneticAlgorithm(self,
                         stop_crit: bool = True,
                         stop_value: float = 1.0)
        '''
        Genetic algorithm
        stop_crit -- True if the algorithm should be stopped when a certain value is reached
        stop_value -- the value for which the algorithm stops, when stop_crit = True
        '''

        # Initial population
        population = [self.random_satellite() for _ in range(self.constellation_size)]
        fitness = [self.evaluation(constellation) for constellation in population]
        print("Initial:", "min = ", np.min(fitness), ", max = ", np.max(fitness), ", average = ", np.mean(fitness))

        # Main loop. Stop criterion: number of iterations
        for num_it in range(self.numIterations):
            if stop_crit and (
                    (self.minimize and np.min(fitness) == stop_value) or (not self.minimize and np.max(fitness) == stop_value)):
                break

            offspring = []

            for kk in range(self.numParents // 2):
                # Parent selection
                parent_1 = population[self.select_parents(fitness)]
                parent_2 = population[self.select_parents(fitness)]

                # Crossover
                child_1, child_2 = self.cross(parent_1, parent_2)

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
            print("Iteration ", num_it, ": min = ", np.min(fitness), ", max = ", np.max(fitness), ", average = ",
                  np.mean(fitness))

        arg = np.argmin(fitness) if self.minimize else np.argmax(fitness)

        return population[arg]