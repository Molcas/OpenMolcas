# Collection of crossover functions for orderings (permutations).
import random
import numpy as np

def order_co(parent1, parent2):
    gene1 = int(random.random() * len(parent1))
    gene2 = gene1
    while gene2 == gene1:
        gene2 = int(random.random() * len(parent1))

    start_gene = min(gene1, gene2)
    end_gene = max(gene1, gene2)

    section1 = parent1[start_gene:end_gene]
    section2 = [gene for gene in parent2 if gene not in section1]

    child = section2[:start_gene] + section1 + section2[start_gene:]

    return child

def my_co(parent1, parent2):
    start_gene = 0
    end_gene = random.randint(1, len(parent1) - 1)

    section1 = parent1[start_gene:end_gene]
    section2 = [gene for gene in parent2 if gene not in section1]

    child = section1 + section2

    return child

def my_co2(parent1, parent2):
    start_gene = 0
    end_gene = random.randint(1, len(parent1) - 1)

    section1 = parent1[start_gene:end_gene]
    if random.random() < 0.5:
        section1.reverse()
    section2 = [gene for gene in parent2 if gene not in section1]

    child = section1 + section2

    return child

def random_co(parent1, parent2):
    return random.sample(parent1, len(parent1))

def partially_mapped_co(parent1, parent2):
    parent1 = np.array(parent1)
    parent2 = np.array(parent2)
    child = np.full(len(parent1), -1, dtype=int)

    gene1 = int(random.random() * len(parent1))
    gene2 = gene1
    while gene2 == gene1:
        gene2 = int(random.random() * len(parent1))

    start_gene = min(gene1, gene2)
    end_gene = max(gene1, gene2)

    child[start_gene:end_gene] = parent1[start_gene:end_gene]
    for i in range(start_gene, end_gene):
        if parent2[i] not in child:
            j = i
            while j >= start_gene and j < end_gene:
                j = np.where(parent2 == parent1[j])[0][0]
            child[j] = parent2[i]

    j = 0
    for i in range(len(child)):
        if child[i] == -1:
            while parent2[j] in child:
                j += 1
            child[i] = parent2[j]
            j += 1

    return list(child)

def shuffle_cluster(parent1, parent2):
    """
    This is not actually a crossover algorithm, as it uses only one parent.
    Clusters the sequence into random-sized groups, reorders the clusters randomly,
    and potentially reverses each cluster with 50% probability.
    
    Args:
        parent1 (list): The sequence to be clustered and reordered
        parent2 (list): Unused parameter, kept for consistency
        
    Returns:
        list: The modified sequence after clustering operations
    """
    length = len(parent1)
    
    num_clusters = random.randint(2,length)

    breakpoints = [0]
    breakpoints.extend(sorted(random.sample(range(1, length), num_clusters - 1)))
    breakpoints.append(length)

    result = []
    clusters = []

    # Create clusters
    for i in range(len(breakpoints) - 1):
        start, end = breakpoints[i], breakpoints[i + 1]
        cluster = parent1[start:end]
        if random.random() < 0.10: # inverse probability is currently hardcoded
            cluster = cluster[::-1]
        clusters.append(cluster)

    # Swap only two clusters instead of full shuffle
    idx1, idx2 = random.sample(range(len(clusters)), 2)
    clusters[idx1], clusters[idx2] = clusters[idx2], clusters[idx1]

    for cluster in clusters:
        result.extend(cluster)

    return result