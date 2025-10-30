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

# This module is an interface between a pre-calculated dataset of GUGA orderings
# and fitnesses in a pandas dataframe. The dataframe is converted into a dictionary
# to utilize a hash table.
# ** In: Dataframe -> process_df.py -> Out: Dictionary **
#
# Outside of this module, the code should only use fitness dictionaries {ordering: fitness}.
#
# We assume an input dataframe has the following structure:
# - The dataframe has a permutation column. The column comprises of orderings that are tuples of integers.
# - The dataframe has fitness columns consisting of (probably floating point) numbers.
# - Column names may differ and the user is supposed to use the column names as input parameters

def gen_fitness_ht(df, perm_col, fitness_col):
    """
    Return a dictionary (hash table) of a selected fitness metric for each permutation.
    - df: A pandas DataFrame of orderings and their fitness metrics
    - perm_col: The name of the column containing the permutations
    - fitness_col: The name of the column containing the fitness scores
    """
    fitness_dict = {tuple(perm[perm_col]): perm[fitness_col] for _, perm in df.iterrows()}

    return fitness_dict

def gen_pop_fitness_ht(population, fitness_dict):
    """
    From a fitness dictionary, generate a reduced fitness dictionary for a given population.
    This mini dictionary is introduced to generalize the code for dynamic fitness evaluation,
    rather than reading everything from df.
    """
    reduced_fitness = {key: fitness_dict[key] for key in population}

    return reduced_fitness

