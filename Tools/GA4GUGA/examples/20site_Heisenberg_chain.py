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
Perform a GA simulation with the S-Ms mapping for a 20-site Heisenberg chain.
`ref_csf` defines the S-Ms consistent CSF.
Therefore, the S-Ms mappling {site:coupling} is {1:1, 2:2, 3:1, 4:2, ..., 19:1, 20:2}.
"""
import GA_mod.run_GA as ga                                                              
import GA_mod.crossover as co                                                           
from GA_mod.measure_fitness import FitnessFunction                                      
                                                                                 
FILEPATH = "../extra_files/"
fcidump = FILEPATH + "FCIDUMP_20site_Heisenberg"

ref_csf = [1,2,1,2,1,2,1,2,1,2,
           1,2,1,2,1,2,1,2,1,2]

ga.perform_GA(
        FitnessFunction.DIAG_ELEM_SMS_MAPPING,
        num_chroms=100,
        restricted_ordering_len=20,
        elite_size=10,
        mutation_rates=[0.01],
        generations=1000,
        co_function=co.order_co,
        fcidump=fcidump,
        norb=20,
        sms_ref_csf=ref_csf,
        tMinimize=True
        )
