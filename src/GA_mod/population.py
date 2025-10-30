import random
from GA_mod import config
from GA_mod import extend_ordering

class Population:
    def __init__(self, num_chroms, ordering_len, elite_size, tExtendChrom,
                 sms_ref_csf=None, sms_ref_ordering=None, restart_filename=None):
        self.num_chroms = num_chroms
        self.ordering_len = ordering_len
        self.elite_size = elite_size
        self.tExtendChrom = tExtendChrom
        # sms_ref_csf assumes the stepvector notation
        self.sms_mapping_enabled = sms_ref_csf is not None and sms_ref_ordering is not None
        if self.sms_mapping_enabled:
            self.sms_mapping_dict = {k:v for k,v in zip(sms_ref_ordering, sms_ref_csf)}
        else:
            self.sms_mapping_dict = None
        if restart_filename is not None:
            self.current_pop = self.read_population(restart_filename)
        else:
            self.current_pop = self.gen_pop_random()

    def gen_pop_random(self):
        """
        Return a list of num_chroms random orderings.
        """
        population = []
        while len(population) < self.num_chroms:
            chrom = tuple(random.sample(range(1, self.ordering_len + 1),
                                        self.ordering_len))
            if ((self.sms_mapping_enabled and is_csf_valid(chrom, self.sms_mapping_dict, self.tExtendChrom))
                or not self.sms_mapping_enabled):
                population.append(chrom)

        return population

    def mutate(self, mutation_rate):
        # mutate only non-elite chromosomes
        for i in range(self.elite_size, self.num_chroms):
            for gene1 in range(0, self.ordering_len):
                if random.random() < mutation_rate:
                    mutated_chrom = list(self.current_pop[i])
                    gene2 = random.sample(range(0, self.ordering_len), 1)[0]
                    mutated_chrom[gene1], mutated_chrom[gene2]\
                  = mutated_chrom[gene2], mutated_chrom[gene1]

                    # prevent invalid CSF
                    if self.sms_mapping_enabled:
                        while not is_csf_valid(mutated_chrom, self.sms_mapping_dict, self.tExtendChrom):
                            mutated_chrom = list(self.current_pop[i])
                            gene2 = random.sample(range(0, self.ordering_len), 1)[0]
                            mutated_chrom[gene1], mutated_chrom[gene2]\
                          = mutated_chrom[gene2], mutated_chrom[gene1]
                    self.current_pop[i] = tuple(mutated_chrom)

    def get_elite_chromosomes(self, reduced_fitness):
        """
        From the current population, get a list of chromosomes with the highest
        fitness.
        """
        elites = []
        sorted_chromosomes = sorted(self.current_pop,
                key=lambda x: reduced_fitness[x], reverse=True)

        self.current_pop = sorted_chromosomes

        for i in range(len(self.current_pop)):
            if self.current_pop[i] not in elites:
                elites.append(self.current_pop[i])
                if len(elites) == self.elite_size:
                    break
        return elites
    
    def gen_mating_pool(self, reduced_fitness, sampling_function):
        """
        Generate a pool of chromosomes for crossover.
        """
        pool = []
        for _ in range(0, self.num_chroms - self.elite_size):
            pool.append(sampling_function(self.current_pop, reduced_fitness))
        return pool
    
    def next_generation(self, reduced_fitness, crossover_function,
                        sampling_function, mutation_rate):
        elites = self.get_elite_chromosomes(reduced_fitness)
        pool = self.gen_mating_pool(reduced_fitness, sampling_function)
        offspring = self.perform_crossover(pool, crossover_function)
        self.current_pop = elites + offspring
        self.mutate(mutation_rate)

    def perform_crossover(self, pool, crossover_function):
        """
        Perform crossover on the given population (list of tuples).
        """
        random.shuffle(pool)

        pool = [list(individual) for individual in pool]

        offsprings = []

        for i in range(0, len(pool)):
            crossover_counter = 0
            max_attempts = 10
            while crossover_counter < max_attempts:
                offspring = crossover_function(pool[i], pool[len(pool)-i-1])
                if ((self.sms_mapping_enabled
                     and is_csf_valid(offspring, self.sms_mapping_dict, self.tExtendChrom))
                    or not self.sms_mapping_enabled):
                    offsprings.append(tuple(offspring))
                    break
                crossover_counter += 1
            if crossover_counter == max_attempts:
                offspring = random.sample(range(1, self.ordering_len + 1),
                                          self.ordering_len)
                while not is_csf_valid(offspring, self.sms_mapping_dict, self.tExtendChrom):
                    offspring = random.sample(range(1, self.ordering_len + 1),
                                              self.ordering_len)
                offsprings.append(tuple(offspring))

        return offsprings

    def read_population(self, filename):
        """
        Initialize population from a file.
        
        File format:
        - Lines beginning with '#' are ignored as comments
        - Each line contains a chromosome in tuple format e.g. (1,3,2,5,...)
        - The last column contain the fitness score (which is ignored)
        """
        chromosomes = []
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') or not line:
                    continue
                    
                # Extract tuple part from the line
                tuple_part = line.split('(')[1].split(')')[0] if '(' in line else line
                
                # Parse chromosome values, handling various formats
                if ',' in tuple_part:
                    # Handle tuple format (1,2,3,...)
                    values = [int(x.strip()) for x in tuple_part.split(',')]
                else:
                    # Handle space-separated format
                    values = [int(x) for x in tuple_part.split()]

                chromosome = tuple(values)
                    
                # If the last value is significantly different, it might be the fitness score
                if len(chromosome) != self.ordering_len:
                    raise ValueError(f"Chromosome length {len(chromosome)} "
                                     f"doesn't match expected length {self.ordering_len}")
                if self.sms_mapping_enabled and not is_csf_valid(chromosome, self.sms_mapping_dict, self.tExtendChrom):
                    raise ValueError(f"Chromosome {chromosome} is not CSF compatible")
                    
                chromosomes.append(chromosome)
        
        # Sanity checks
        if len(chromosomes) != self.num_chroms:
            raise ValueError(f"Number of chromosomes in file ({len(chromosomes)}) "
                             f"doesn't match expected number ({self.num_chroms})")
        
        return chromosomes

#-------------------------------------------------------------------------------

def is_csf_valid(ordering, sms_mapping_dict, tExtendChrom):
    if tExtendChrom:
        new_ordering = extend_ordering.extend_ordering(
            ordering, 
            config.on_site_permutation, 
            config.num_prefix, 
            config.num_suffix
        )
    else:
        new_ordering = ordering

    cumul_spin = 0
    for i in new_ordering:
        if sms_mapping_dict[i] == 1:
            cumul_spin += 1
        elif sms_mapping_dict[i] == 2:
            cumul_spin -= 1

        if cumul_spin < 0:
            return False
    return True
