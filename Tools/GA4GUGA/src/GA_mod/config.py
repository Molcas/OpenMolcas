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
Global configuration for GA_mod package.

These variables are set by perform_GA and can be accessed by other modules.
"""

# Ordering extension parameters
on_site_permutation = (1,)
num_prefix = 0
num_suffix = 0

def set_ordering_params(on_site_perm, prefix, suffix):
    """
    Set global ordering parameters.
    
    Args:
        on_site_perm: Tuple of integers representing the on-site pattern to repeat
        prefix: Number of orbitals to prepend
        suffix: Number of orbitals to append
    """
    global on_site_permutation, num_prefix, num_suffix
    on_site_permutation = on_site_perm
    num_prefix = prefix
    num_suffix = suffix
