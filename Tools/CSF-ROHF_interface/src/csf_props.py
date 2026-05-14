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
# Copyright (C) 2026, Maru Song                                        *
#***********************************************************************

"""CSF-specific helper functions for MOLCAS validation."""

from __future__ import annotations

import re
from typing import Iterable


def calculate_csf_active_electrons(csf_stepvec: Iterable[int]) -> int:
    """Calculate the active electron count implied by a CSF step vector."""
    electron_map = {
        0: 0,
        1: 1,
        2: 1,
        3: 2,
    }

    active_electrons = 0
    for step in csf_stepvec:
        try:
            step_value = int(step)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid CSF step vector element: {step!r}") from exc

        if step_value not in electron_map:
            raise ValueError(f"Unsupported CSF step vector element: {step_value}")

        active_electrons += electron_map[step_value]

    return active_electrons


def calculate_csf_spin_twice(csf_stepvec: Iterable[int]) -> int:
    """Calculate 2S implied by a CSF step vector.

    The cumulative value must never become negative while scanning left to right.
    """
    spin_map = {
        0: 0,
        1: 1,
        2: -1,
        3: 0,
    }

    spin_twice = 0
    for step in csf_stepvec:
        try:
            step_value = int(step)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"Invalid CSF step vector element: {step!r}") from exc

        if step_value not in spin_map:
            raise ValueError(f"Unsupported CSF step vector element: {step_value}")

        spin_twice += spin_map[step_value]
        if spin_twice < 0:
            raise ValueError(
                f"Invalid CSF step vector: cumulative 2S became negative ({spin_twice})"
            )

    return spin_twice


def extract_active_orbitals(log_content: str) -> int:
    """Extract the active orbital count from MOLCAS log content."""
    matches = re.findall(r"Number of active orbitals\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'Number of active orbitals' in MOLCAS log")

    return int(matches[-1])


def extract_active_shell_electrons(log_content: str) -> int:
    """Extract the active-shell electron count from MOLCAS log content."""
    matches = re.findall(r"Number of electrons in active shells\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'Number of electrons in active shells' in MOLCAS log")

    return int(matches[-1])


def extract_state_symmetry(log_content: str) -> int:
    """Extract the state symmetry / multiplicity value from MOLCAS log content."""
    matches = re.findall(r"State symmetry\s+(\d+)", log_content)
    if not matches:
        raise ValueError("Could not find 'State symmetry' in MOLCAS log")

    return int(matches[-1])


def validate_csf_against_molcas(
    log_content: str,
    csf_orbital_count: int,
    csf_active_electrons: int,
    csf_spin_twice: int,
) -> None:
    """Validate precomputed CSF properties against the first-RASSCF MOLCAS log information."""

    active_orbitals = extract_active_orbitals(log_content)
    if csf_orbital_count != active_orbitals:
        raise ValueError(
            f"CSF orbital count ({csf_orbital_count}) does not match the number of active orbitals ({active_orbitals})"
        )

    active_shell_electrons = extract_active_shell_electrons(log_content)
    if csf_active_electrons != active_shell_electrons:
        raise ValueError(
            f"CSF electron count ({csf_active_electrons}) does not match the number of electrons in active shells "
            f"({active_shell_electrons})"
        )

    state_symmetry = extract_state_symmetry(log_content)
    csf_multiplicity = csf_spin_twice + 1
    if csf_multiplicity != state_symmetry:
        raise ValueError(
            f"CSF multiplicity ({csf_multiplicity}) does not match MOLCAS spin symmetry ({state_symmetry})"
        )
