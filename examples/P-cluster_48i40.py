"""
Perform a GA simulation for CAS(48,40) PN-cluster.
The fitness function maximizes the gap between the highest and the lowest
CSF energies over the 14 collinear CSFs (`collinear_csfs`).

Since we restrict the permutation space to only permutations among the 8 Fe
sites, `restricted_ordering_len` is set to 8, and `on_site_permutation`
is set to (1,2,3,4,5) meaning there is no permutation within each Fe site.
"""
import GA_mod.run_GA as ga                                                              
import GA_mod.crossover as co                                                           
from GA_mod.measure_fitness import FitnessFunction                                      

fcidump = '../extra_files/FCIDUMP_PN_48i40'
collinear_csfs = [
[3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2],
[3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,1,1,1,1,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2],
[3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,1,1,1,1,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2,3,2,2,2,2]
        ]

ga.perform_GA(
        FitnessFunction.MIN_MAX_DIFF,
        num_chroms=20,
        elite_size=2,
        mutation_rates=[0.02],
        generations=1000,
        restricted_ordering_len=8,
        on_site_permutation=(1,2,3,4,5),
        num_suffix=0,
        co_function=co.order_co,
        fcidump=fcidump,
        norb=40,
        csf_list=collinear_csfs,
        )
