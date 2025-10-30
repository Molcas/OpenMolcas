import sys
import os
from GA_mod import population as pop
from GA_mod import sampling
from GA_mod import measure_fitness
from GA_mod import crossover as co
from GA_mod import extend_ordering
from GA_mod import config
from FCIDUMP_tools import IntegralClass
import subprocess

def perform_GA(fitness_function, num_chroms, elite_size, mutation_rates,
               restricted_ordering_len, generations, co_function, fcidump, norb,
               cluster_period=5, stagnation_limit=100,
               on_site_permutation=(1,), num_prefix=0, num_suffix=0, 
               restart_filename=None, **kwargs):
    """
    Main function to perform Genetic Algorithm optimization.
    """


    expected_norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)
    if expected_norb != norb:
        raise ValueError(
            "Inconsistent norb: num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation) "
            f"= {expected_norb}, but norb = {norb}"
        )
    if norb == restricted_ordering_len:
        tExtendChrom = False
    else:
        tExtendChrom = True

    # Set global ordering parameters for use in other modules
    # TODO: set other variables as global as well.
    config.set_ordering_params(on_site_permutation, num_prefix, num_suffix)

    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD']).decode('ascii').strip()

    pop_filename = kwargs.get('pop_file_name', 'current_pop.log')

    checkpoint_trigger = kwargs.get('checkpoint_trigger', 'WRITE_CHECKPOINT')
    checkpoint_prefix = kwargs.get('checkpoint_prefix', 'FCIDUMP_checkpoint')

    sms_ref_csf = kwargs.get('sms_ref_csf', None)

    print("Genetic Algoritm simulation started", file=sys.stdout)
    print(f"GIT hash: {git_hash}", file=sys.stdout)
    print("", file=sys.stdout)
    if restart_filename is not None:
        print(f"Restarting from population file: {restart_filename}\n", file=sys.stdout)
    print(f"- Number of chromosomes: {num_chroms}", file=sys.stdout)
    print(f"- Elite size: {elite_size}", file=sys.stdout)
    print(f"- Mutation rates: {mutation_rates}", file=sys.stdout)
    print(f"- Generations: {generations}", file=sys.stdout)
    print(f"- Crossover function: {co_function.__name__}", file=sys.stdout)
    print(f"- Fitness function: {fitness_function}", file=sys.stdout)
    print(f"- Clustering period: {cluster_period}", file=sys.stdout)
    print(f"- Stagnation limit: {stagnation_limit}", file=sys.stdout)
    if sms_ref_csf is not None:
        print(f"- S-Ms mapping reference CSF: {sms_ref_csf}", file=sys.stdout)
        sms_ref_ordering = tuple(range(1, norb + 1))
    else:
        sms_ref_ordering = None
    print("", file=sys.stdout)
    print(f"Checkpoint trigger: Create file '{checkpoint_trigger}' to write current best ordering", file=sys.stdout)
    print("", file=sys.stdout)

    POPClass = pop.Population(num_chroms, restricted_ordering_len, elite_size,
                               tExtendChrom, sms_ref_csf, sms_ref_ordering,
                               restart_filename=restart_filename)

    FCIDUMPClass = IntegralClass.FCIDUMPReader(fcidump)

    # Track best fitness and stagnation
    stagnation_counter = 0
    mutation_rate_index = 0
    current_mutation_rate = mutation_rates[mutation_rate_index]

    # 0th generation
    reduced_fitness_dict = \
        measure_fitness.calculate_fitness(fitness_function, POPClass, FCIDUMPClass,
                                          norb, tExtendChrom=tExtendChrom,
                                          **kwargs)
    bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)

    best_fitness = reduced_fitness_dict[bestchrom]
    print("# Generation  Ordering  Fitness", file=sys.stdout)
    extended_bestchrom = extend_ordering.extend_ordering(bestchrom, 
        on_site_permutation, num_prefix, num_suffix)
    print(f"0  {extended_bestchrom}  {best_fitness}", file=sys.stdout)

    # Subsequent generations
    for i in range(1, generations + 1):
        best_fitness_prev = best_fitness
        if i % cluster_period == 0:
            _co_function = co.shuffle_cluster
        else:
            _co_function = co_function
        POPClass.next_generation(reduced_fitness_dict, _co_function,
                                  sampling.roullette_wheel_sampling,
                                  current_mutation_rate)

        reduced_fitness_dict = \
            measure_fitness.calculate_fitness(fitness_function, POPClass,
                                              FCIDUMPClass, norb,
                                              tExtendChrom=tExtendChrom,
                                              **kwargs)
        bestchrom = max(reduced_fitness_dict, key=reduced_fitness_dict.get)
        best_fitness = reduced_fitness_dict[bestchrom]

        if abs(best_fitness - best_fitness_prev) < 1e-6:
            stagnation_counter += 1
            if stagnation_counter == stagnation_limit:
                mutation_rate_index = (mutation_rate_index + 1) % len(mutation_rates)
                current_mutation_rate = mutation_rates[mutation_rate_index]
                print(f"# Stagnation detected. Change mutation rate to {current_mutation_rate}", file=sys.stdout)
                stagnation_counter = 0
        else:
            stagnation_counter = 0

        extended_bestchrom = extend_ordering.extend_ordering(bestchrom, 
            on_site_permutation, num_prefix, num_suffix)
        print(f"{i}  {extended_bestchrom}  {best_fitness}", file=sys.stdout)

        # Check for checkpoint trigger file
        if os.path.exists(checkpoint_trigger):
            checkpoint_filename = f"{checkpoint_prefix}_gen{i}"
            print(f"\n# Checkpoint trigger detected at generation {i}", file=sys.stdout)
            print(f"# Writing FCIDUMP with current best ordering to '{checkpoint_filename}'", file=sys.stdout)
            FCIDUMPClass.dump_integrals(checkpoint_filename, extended_bestchrom)
            print(f"# Checkpoint written. Removing trigger file and continuing...\n", file=sys.stdout)
            try:
                os.remove(checkpoint_trigger)
            except OSError as e:
                print(f"# Warning: Could not remove trigger file: {e}", file=sys.stdout)

        with open(pop_filename, 'w') as log_file:
            log_file.write(f"# Chromosomes in the {i}th generation and their fitnesses\n")
            for chrom in POPClass.current_pop:
                extended_chrom = extend_ordering.extend_ordering(chrom, 
                    on_site_permutation, num_prefix, num_suffix)
                log_file.write(f"{extended_chrom} {reduced_fitness_dict[chrom]}\n")

    # Generate an FCIDUMP file with the best ordering.
    print("\n\nGenetic Algorithm simulation completed.", file=sys.stdout)
    print(f"Best ordering found: {extended_bestchrom}", file=sys.stdout)
    print("Writing reordered FCIDUMP to 'FCIDUMP_bestordering'", file=sys.stdout)
    FCIDUMPClass.dump_integrals('FCIDUMP_bestordering', extended_bestchrom)
