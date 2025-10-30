"""
Perform a GA simulation with the S-Ms mapping for a 20-site Heisenberg chain.
`ref_csf` defines the S-Ms consistent CSF.
Therefore, the S-Ms mappling {site:coupling} is {1:1, 2:2, 3:1, 4:2, ..., 19:1, 20:2}.

`FAST_DIAG_MIN_OSONL` "minimizes" the (S-MS consistent) CSF energy using only
ordering dependent terms for open-shell orbitals (1 and 2).
This fitness function only works with Heisenberg Hamiltonians, as it minimizes
the CSF energy and assumes CSFs only consist of open-shell orbitals.

The code first creats a Heisenberg FCIDUMP file and then runs the GA simulation with
the FCIDUMP file.
"""
import GA_mod.run_GA as ga                                                              
import GA_mod.crossover as co                                                           
from GA_mod.measure_fitness import FitnessFunction                                      
import FCIDUMP_tools.xyz2heisenberg as x2h
                                                                                 
FILEPATH = "../extra_files/"
xyzfile = FILEPATH + "20-site-chain.xyz"
fcidump = 'FCIDUMP_20site_Heisenberg_chain_nn'

ref_csf = [1,2,1,2,1,2,1,2,1,2,
           1,2,1,2,1,2,1,2,1,2]

# Create a Heisenberg FCIDUMP file where J = 1/r using the geometry in 
# `xyzfile` only considering 0.9 =< r =< 1.1.
x2h.gen_heisenberg_fcidump(xyzfile, fcidump, nel=20, norb=20, ms=0, distance_range=[0.9,1.1])

ga.perform_GA(
        FitnessFunction.FAST_DIAG_MIN_OSONLY,
        num_chroms=200,
        restricted_ordering_len=20,
        elite_size=20,
        mutation_rates=[0.01],
        generations=10000,
        cluster_period=2,
        co_function=co.order_co,
        fcidump=fcidump,
        norb=20,
        sms_ref_csf=ref_csf
        )

