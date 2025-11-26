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
Wrapper to measure fitness scores.
The larger the score, the better the fitness.
The sign of the fitness doesn't matter, as it is rescaled in sampling.py
But when targetting a minimum, -1 has to be multipled at the end to give the
minimum a larger score.
"""
from enum import Enum, auto
from GA_mod import extend_ordering
from GA_mod import GUGA_diag
from GA_mod import config
from GA_mod import auxiliary_functions as af
import numpy as np

class FitnessFunction(Enum):
    MIN_MAX_DIFF = auto()
    MAX_DIAG_ELEM = auto()
    MIN_DIAG_ELEM = auto()
    DIAG_ELEM_SMS_MAPPING = auto()
    FAST_DIAG_MIN_OSONLY = auto()

def _diag_elem_sms_mapping(population, extended_pop, sms_mapping_dict, FCIDUMPClass, norb, tMinimize):
    """Internal function that calculates fitness based on reference diagonal elements."""
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        CSF_ref = [sms_mapping_dict[orb] for orb in extended_chrom]

        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)

        diagelem = GUGAClass.calc_diag_elem(CSF_ref, add_core=True)

        if tMinimize:
            reduced_fitness[chrom] = abs(diagelem)
        else:
            reduced_fitness[chrom] = diagelem
    return reduced_fitness

def _min_max_diff_fitness(population, extended_pop, FCIDUMPClass, norb, csf_list):
    """Internal function that calculates fitness based on min-max difference."""
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        min_val = float('inf')
        max_val = float('-inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            min_val = min(min_val, diagelem)
            max_val = max(max_val, diagelem)
        reduced_fitness[chrom] = max_val - min_val
    return reduced_fitness

def _max_diagelem(population, extended_pop, FCIDUMPClass, norb, csf_list):
    """
    Internal function that calculates fitness based the max diagonal element
    using csf_list.
    """
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        max_val = float('-inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            max_val = max(max_val, diagelem)
        reduced_fitness[chrom] = max_val
    return reduced_fitness

def _min_diagelem(population, extended_pop, FCIDUMPClass, norb, csf_list):
    """
    Internal function that calculates fitness based the min diagonal element
    using csf_list.
    """
    reduced_fitness = {}
    for chrom, extended_chrom  in zip(population, extended_pop):
        FCIDUMPClass.permute_integrals(extended_chrom, t_passive=True)
        GUGAClass = GUGA_diag.DiagElement(norb, FCIDUMPClass)
        min_val = float('inf')
        for csf_stepvec in csf_list:
            diagelem = GUGAClass.calc_diag_elem(csf_stepvec, add_core=True)
            min_val = min(min_val, diagelem)
        reduced_fitness[chrom] = min_val * -1
    return reduced_fitness

def _fast_diag_min_osonly(population, extended_pop, J, sms_mapping_dict):
    reduced_fitness = {}
    J = np.array(J)
    for chrom, extended_chrom  in zip(population, extended_pop):
        sorting_arr = np.array(extended_chrom) - 1
        J_reordered = J[np.ix_(sorting_arr, sorting_arr)]

        csf_stepvec = [sms_mapping_dict[i] for i in extended_chrom]
        X = af.X_matrix_openshell_only(csf_stepvec)

        reduced_fitness[chrom] = np.sum(J_reordered * X) * -1
    return reduced_fitness

#------------------------------------------------------------------------------#

def calculate_fitness(method: FitnessFunction, POPClass, FCIDUMPClass, norb, 
                      tExtendChrom, **kwargs):
    """
    Wrapper function to calculate fitness using specified method.

    Args:
        method (FitnessFunction): The method to use for fitness score evaluation
        POPClass: Population class instance
        FCIDUMPClass: FCIDUMP class instance
        norb: Number of orbitals
        tExtendChrom: Whether to extend chromosomes using ordering parameters
        **kwargs: Optional arguments depending on method:
            - For MIN_MAX_DIFF: 
                csf_list (list): List of CSFs to measure
            - For future methods:
                Add parameters as needed

    Returns:
        dict: Dictionary mapping chromosomes to their fitness values

    Raises:
        ValueError: If required parameters are missing for chosen method
    """
    sms_mapping_dict = POPClass.sms_mapping_dict
    csf_list = kwargs.get('csf_list', None)
    tMinimize = kwargs.get('tMinimize', False)

    if tExtendChrom:
        extended_pop = [extend_ordering.extend_ordering(ordering, 
            config.on_site_permutation, config.num_prefix, config.num_suffix) 
            for ordering in POPClass.current_pop]
    else:
        extended_pop = POPClass.current_pop

    if method == FitnessFunction.MIN_MAX_DIFF:
        if csf_list is None:
            raise ValueError("csf_list is required for MIN_MAX_DIFF method")
        fitness_ht = _min_max_diff_fitness(POPClass.current_pop, extended_pop, FCIDUMPClass, norb, csf_list)
    elif method == FitnessFunction.DIAG_ELEM_SMS_MAPPING:
        if sms_mapping_dict is None:
            raise ValueError("sms_mapping_dict is required for FAST_DIAG_MIN_OSONLY method")
        fitness_ht = _diag_elem_sms_mapping(POPClass.current_pop, extended_pop, sms_mapping_dict, FCIDUMPClass, norb, tMinimize)
    elif method == FitnessFunction.MAX_DIAG_ELEM:
        fitness_ht = _max_diagelem(POPClass.current_pop, extended_pop, FCIDUMPClass, norb, csf_list)
    elif method == FitnessFunction.MIN_DIAG_ELEM:
        fitness_ht = _min_diagelem(POPClass.current_pop, extended_pop, FCIDUMPClass, norb, csf_list)
    elif method == FitnessFunction.FAST_DIAG_MIN_OSONLY:
        if sms_mapping_dict is None:
            raise ValueError("sms_mapping_dict is required for FAST_DIAG_MIN_OSONLY method")
        J = af.J_mat_from_fcidump(FCIDUMPClass, norb)
        fitness_ht = _fast_diag_min_osonly(POPClass.current_pop, extended_pop, J, sms_mapping_dict)
    else:
        raise ValueError(f"Unknown fitness method: {method}")

    return fitness_ht
