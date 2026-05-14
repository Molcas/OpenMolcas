#!/usr/bin/env python3
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
Main executable script for running GA4GUGA genetic algorithm optimization.
This script provides a command-line interface for the GA4GUGA package.
"""

import sys
import os
import argparse
import ast

# Add the src directory to Python path to import GA_mod
# GA_mod is located in src/GA_mod relative to this script
script_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(script_dir, 'src')
if src_dir not in sys.path:
    sys.path.insert(0, src_dir)

import GA_mod.run_GA as ga
import GA_mod.crossover as co
from GA_mod.measure_fitness import FitnessFunction


def parse_list_arg(arg_str):
    """Parse a string representation of a list into an actual list."""
    try:
        result = ast.literal_eval(arg_str)
        if isinstance(result, list):
            return result
        else:
            raise ValueError("Not a list")
    except:
        raise argparse.ArgumentTypeError(f"Invalid list format: {arg_str}")


def parse_tuple_arg(arg_str):
    """Parse a string representation of a tuple into an actual tuple."""
    try:
        result = ast.literal_eval(arg_str)
        if isinstance(result, tuple):
            return result
        else:
            raise ValueError("Not a tuple")
    except:
        raise argparse.ArgumentTypeError(f"Invalid tuple format: {arg_str}")


def parse_csf_list_arg(arg_str):
    """Parse a string representation of a nested list (CSF list) into an actual list of lists."""
    try:
        result = ast.literal_eval(arg_str)
        if isinstance(result, list) and all(isinstance(item, list) for item in result):
            return result
        else:
            raise ValueError("Not a list of lists")
    except:
        raise argparse.ArgumentTypeError(f"Invalid CSF list format: {arg_str}")


def main():
    parser = argparse.ArgumentParser(
        description='Run genetic algorithm optimization for GUGA orbital ordering.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # 20-site Heisenberg chain optimization with S-Ms mapping
    ga_driver.py \
      --fcidump $path_to_extrafiles/FCIDUMP_20site_Heisenberg \
      --norb 20 \
      --fitness DIAG_ELEM_SMS_MAPPING \
      --sms-ref-csf "[1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2]" \
      --num-chroms 100 \
      --elite-size 10 \
      --mutation-rates "[0.01]" \
      --generations 1000 \
      --restricted-ordering-len 20 \
      --minimize
        """
    )

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('--fcidump', type=str, required=True,
                          help='Path to the FCIDUMP file')
    required.add_argument('--norb', type=int, required=True,
                          help='Total number of orbitals. Must satisfy: '
                               'norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)')
    required.add_argument('--fitness', type=str, required=True,
                          choices=['MIN_MAX_DIFF', 'MAX_DIAG_ELEM', 'MIN_DIAG_ELEM', 
                                   'DIAG_ELEM_SMS_MAPPING', 'FAST_DIAG_MIN_OSONLY'],
                          help='Fitness function: MIN_MAX_DIFF (maximize diff between highest/lowest CSF energies), '
                               'MAX_DIAG_ELEM (maximize highest CSF energy), MIN_DIAG_ELEM (minimize lowest CSF energy), '
                               'DIAG_ELEM_SMS_MAPPING (max/min S-Ms consistent CSF energy), '
                               'FAST_DIAG_MIN_OSONLY (minimize S-Ms consistent CSF energy, evaluate ordering-dependent terms only)')
    required.add_argument('--num-chroms', type=int, required=True,
                          help='Population size (number of chromosomes)')
    required.add_argument('--elite-size', type=int, required=True,
                          help='Number of elite chromosomes retained each generation')
    required.add_argument('--mutation-rates', type=parse_list_arg, required=True,
                          help='List of mutation rates, e.g., "[0.01]" or "[0.01, 0.02]". '
                               'If multiple rates provided, algorithm cycles through them when progress stalls')
    required.add_argument('--generations', type=int, required=True,
                          help='Number of generations to run')
    required.add_argument('--restricted-ordering-len', type=int, required=True,
                          help='Length of the effective (restricted) ordering that GA optimizes')

    # Optional arguments
    parser.add_argument('--crossover', type=str, default='order_co',
                        choices=['order_co'],
                        help='Crossover operator for GA simulation. Currently only order crossover is implemented (default: order_co)')
    parser.add_argument('--cluster-period', type=int, default=5,
                        help='Period of cluster shuffling replacing the crossover step (default: 5)')
    parser.add_argument('--stagnation-limit', type=int, default=100,
                        help='If best fitness does not evolve for n generations, '
                             'mutation rate updates to next value in mutation_rates (default: 100)')
    parser.add_argument('--on-site-permutation', type=parse_tuple_arg, default='(1,)',
                        help='Pattern applied to each gene when expanding restricted ordering to full ordering. '
                             'E.g., "(1,2,3,4,5)" for on-site permutation (default: "(1,)")')
    parser.add_argument('--num-prefix', type=int, default=0,
                        help='Number of leading orbitals that remain fixed (not optimized) (default: 0)')
    parser.add_argument('--num-suffix', type=int, default=0,
                        help='Number of trailing orbitals that remain fixed (not optimized) (default: 0)')
    parser.add_argument('--restart-file', type=str, default=None,
                        help='Population filename. If specified, GA will restart from the population file')
    parser.add_argument('--pop-file-name', type=str, default='current_pop.log',
                        help='File to store latest population (set of orderings) during GA run (default: current_pop.log)')
    parser.add_argument('--checkpoint-trigger', type=str, default='WRITE_CHECKPOINT',
                        help='File name whose presence triggers writing FCIDUMP with current best ordering (default: WRITE_CHECKPOINT)')
    parser.add_argument('--checkpoint-prefix', type=str, default='FCIDUMP_checkpoint',
                        help='Prefix for checkpoint FCIDUMP files (default: FCIDUMP_checkpoint)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random number seed for reproducibility (default: None)')

    # Fitness-specific arguments
    fitness_specific = parser.add_argument_group('fitness-specific arguments')
    fitness_specific.add_argument('--sms-ref-csf', type=parse_list_arg, default=None,
                        help='(For DIAG_ELEM_SMS_MAPPING/FAST_DIAG_MIN_OSONLY) '
                             'Set S-Ms consistent CSF in natural ordering (step-vector notation: '
                             '0=empty, 1=spin-up, 2=spin-down, 3=doubly-occupied). '
                             'E.g., "[1,2,1,2]" maps {orbital:coupling} = {1:1, 2:2, 3:1, 4:2}')
    fitness_specific.add_argument('--minimize', action='store_true',
                        help='(For DIAG_ELEM_SMS_MAPPING only) Minimize fitness score instead of maximize')
    fitness_specific.add_argument('--csf-list', type=parse_csf_list_arg, default=None,
                        help='(For MIN_MAX_DIFF/MAX_DIAG_ELEM/MIN_DIAG_ELEM) '
                             'List of CSFs for fitness evaluation (step-vector notation). '
                             'E.g., "[[1,2,1,2],[2,1,2,1]]"')

    args = parser.parse_args()

    expected_norb = args.num_prefix + args.num_suffix + args.restricted_ordering_len * len(args.on_site_permutation)
    if expected_norb != args.norb:
        parser.error(f"Inconsistent norb: num_prefix ({args.num_prefix}) + num_suffix ({args.num_suffix}) + "
                     f"restricted_ordering_len ({args.restricted_ordering_len}) * "
                     f"len(on_site_permutation) ({len(args.on_site_permutation)}) = {expected_norb}, "
                     f"but norb = {args.norb}")

    fitness_func = getattr(FitnessFunction, args.fitness)

    co_func = getattr(co, args.crossover)

    kwargs = {
        'pop_file_name': args.pop_file_name,
        'checkpoint_trigger': args.checkpoint_trigger,
        'checkpoint_prefix': args.checkpoint_prefix,
    }

    if args.sms_ref_csf is not None:
        kwargs['sms_ref_csf'] = args.sms_ref_csf

    if args.minimize:
        kwargs['tMinimize'] = True

    if args.csf_list is not None:
        kwargs['csf_list'] = args.csf_list

    # Run GA
    try:
        ga.perform_GA(
            fitness_function=fitness_func,
            num_chroms=args.num_chroms,
            elite_size=args.elite_size,
            mutation_rates=args.mutation_rates,
            restricted_ordering_len=args.restricted_ordering_len,
            generations=args.generations,
            co_function=co_func,
            fcidump=args.fcidump,
            norb=args.norb,
            cluster_period=args.cluster_period,
            stagnation_limit=args.stagnation_limit,
            on_site_permutation=args.on_site_permutation,
            num_prefix=args.num_prefix,
            num_suffix=args.num_suffix,
            restart_filename=args.restart_file,
            seed=args.seed,
            **kwargs
        )
    except Exception as e:
        print(f"\nError running GA: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
