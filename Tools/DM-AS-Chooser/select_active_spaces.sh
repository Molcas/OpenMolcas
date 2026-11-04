#! /bin/sh

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
# Copyright (C) 2026, Benjamin W. Kaufold                              *
#               2026, Meredith E. Wanta                                *
#               2026, Sijia S. Dong                                    *
#***********************************************************************


#SBATCH --job-name=active_space_chooser
#SBATCH --output=bash_logs/active_space_chooser-%j.log
#SBATCH --error=bash_logs/active_space_chooser-%j.err
#SBATCH --time=30:00
#SBATCH --partition=short,donglab
#SBATCH -C cascadelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --requeue

# --methods: the protocols to use for selecting the active space.
#            format: --methods method1 method2 methodN
#            currently implemented: gdm-as edm-as d2dm-as dir-gdm-as dir-edm-as dir-d2dm-as
# --data-dir: the absolute file path the the directory containing all candidate active spaces and reference calculations.
# --output-dir: the absolute file path to the desired output directory
# -m/--mr-keywords: the active space candidate calculation keywords (e.g. casci, casscf, etc.)
#                   format: -m casci casscf
# -r/--ref-keywords: the reference TD_DFT calculation keywords (e.g. pbe, blyp, etc.)
#                   format: -r pbe blyp
# --extract-excdata: OPTIONAL; to run the excited-state data extraction workflow. Is only available for
#                    CASSCF, CASPT2, and CAS-PDFT
# --generate-comparison-input: OPTIONAL; to generate input for follow-up comparison calculations for selected active spaces.


python3 main.py \
  --methods gdm-as edm-as d2dm-as dir-gdm-as dir-edm-as dir-d2dm-as \
  --data-dir /projects/donglab/wanta.m/open_molcas_integration/sample_data/acetaldehyde \
  --output-dir /projects/donglab/wanta.m/open_molcas_integration/system_results/acetaldehyde \
  -m casci casscf \
  -r b3lyp blyp camb3lyp lcwhpbe m06 m062x m06hf m06l pbe pbe0 r2scan0 r2scan50 r2scanh \
  --extract-excdata \
  --generate-comparison-input
