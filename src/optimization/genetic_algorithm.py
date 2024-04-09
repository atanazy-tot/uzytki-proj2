class genetic_algorithm:
    def __init__(self,
                 ground_stations,
                 planet,
                 start_time,
                 end_time,
                 time_step,
                 numSatellites: int,
                 popSize: int,
                 numParents: int,
                 pm: float,
                 numIterations: int):

        assert numParents % 2 == 0

        self.ground_stations = ground_stations
        self.planet = planet
        self.start_time = start_time
        self.end_time = end_time
        self.time_step = time_step
        self.numSatellites = numSatellites
        self.popSize = popSize
        self.numParents = numParents
        self.pm = pm
        self.numIterations = numIterations
        self.minimize = True

    def evaluation(self, satellites, ground_stations, planet, start_time, end_time, time_step): # In our case, we want to maximize coverage
        self.minimize = False
        c = CoverageSimulator(satellites, ground_stations, planet)
        result = c.calculate_covarage(start_time, end_time, time_step)
        return result

    def mutation(self, satellites):

        pass

    def cross(self):
        pass

    def select_parents(self, fitness, popSize, k = 2): # tournament selection
        """
        Draws k individuals and returns the best one
        fitness -- fitness matrix for population
        popSize -- population size
        k -- number of participants in the tournament.
        """
        candidates = random.sample(range(popSize), k)
        arg = np.argmin(fitness[candidates]) if self.minimize else np.argmax(fitness[candidates])
        return candidates[arg]

    def survival(self,
                 population: np.ndarray,
                 fitness: np.ndarray,
                 offspring: np.ndarray,
                 fitnessOffspring: np.ndarray,
                 popSize: int,
                 permuSize: int) -> np.ndarray:
        """
        finds survivors
        population -- ,
        fitness -- ,
        offspring -- ,
        fitnessOffspring -- ,
        popSize -- population size,
        permuSize -- permutation size
        """
        currentGeneration = np.vstack((population, offspring))
        fitnessGeneration = np.concatenate((fitness, fitnessOffspring))

        order = np.argsort(fitnessGeneration)

        survivors = np.zeros([popSize, permuSize + 1])
        survivors[:, :-1] = currentGeneration[order[:popSize]]
        survivors[:, -1] = fitnessGeneration[order[:popSize]]

        return survivors

    def geneticAlgorithm(self,
                         stop_crit: bool = True,
                         stop_value: float = 1.0) -> np.ndarray:
        '''
        Genetic algorithm
        stop_crit -- True if the algorithm should be stopped when a certain value is reached
        stop_value -- the value for which the algorithm stops, when stop_crit = True
        '''

        # Initial population
        population = []
        fitness = np.array([self.evaluation(satellites, self.ground_stations, self.planet, self.start_time, self.end_time, self.time_step) for satellites in population])
        print("Initial:", "min = ", min(fitness), ", max = ", max(fitness), ", average = ", np.mean(fitness))

        # Main loop. Stop criterion: number of iterations
        for num_it in range(numIterations):
            if stop_crit and (
                    (minimize and min(fitness) == stop_value) or (not minimize and max(fitness) == stop_value)):
                break

            offspringMatrix = np.zeros([numParents, permuSize])

            for kk in range(numParents // 2):
                # Parent selection
                candidate1 = select_parents(fitness, popSize, *args, **kwargs)
                candidate2 = select_parents(fitness, popSize, *args, **kwargs)

                # Crossover
                children = np.array(
                    [cross(population[[candidate1, candidate2]], permuSize, copied_parent=0, *args, **kwargs),
                     cross(population[[candidate1, candidate2]], permuSize, copied_parent=1, *args, **kwargs)])

                # Mutation
                children = np.array([mutation(children[0], pm), mutation(children[0], pm)])

                # Save the children
                offspringMatrix[(kk * 2):(kk * 2 + 2), :] = children

            # Evaluation of children
            fitnessChildren = np.array(
                [evaluation(offspringMatrix[individual,], *args, **kwargs) for individual in range(numParents)])

            # Survivor selection
            res = survival(population, fitness, offspringMatrix, fitnessChildren, popSize, permuSize, *args, **kwargs)
            population = res[:, 0:permuSize]
            fitness = res[:, permuSize]

            # Print info
            print("Iteration ", num_it, ": min = ", min(fitness), ", max = ", max(fitness), ", average = ",
                  np.mean(fitness))

        arg = np.argmin(fitness) if minimize else np.argmax(fitness)

        return population[arg]

    pass