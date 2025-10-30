import random
import numpy as np

def rescale_fitness(fitness):
    min_fitness = min(fitness.values())
    max_fitness = max(fitness.values())

    if max_fitness < 0:
        shift = abs(min_fitness + max_fitness)
        rescaled_fitness = {chromosome:
                        (fitness[chromosome] + shift) for chromosome in fitness}
    else:
        rescaled_fitness = fitness.copy()

    return rescaled_fitness

def roullette_wheel_sampling(population, reduced_fitness):
    """
    Return a chromosome sampled with a roullete wheel sampling.

    One may consider returning multiple chromosomes by setting multiple arrows
    on the whell, rather than doing multiple roullete wheel samplings.
    """
    rescaled_fitness = rescale_fitness(reduced_fitness)
    # rescaled_fitness = reduced_fitness
    cum_sum = 0.0
    cum_arr = np.zeros(len(population))
    for i, chromosome in enumerate(population):
        cum_sum += rescaled_fitness[chromosome]
        cum_arr[i] = cum_sum

    random_num = random.uniform(0, cum_arr[-1])
    idx = np.searchsorted(cum_arr, random_num)
    sampled_chromosome = population[idx]

    return sampled_chromosome
