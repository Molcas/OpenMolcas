#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2025, Maru Song                                        *
#***********************************************************************

"""
Restarts from a population file.
To restart, the population filename is specfied in the `restart_filename` argument.
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
        sms_ref_csf=ref_csf,
        restart_filename=FILEPATH + 'NN_Heisenberg_20chain.pop'
        )

