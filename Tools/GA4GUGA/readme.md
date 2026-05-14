# GA4GUGA

A genetic algorithm (GA) for finding optimal orbital orderings that yield
compact wave function representations in GUGA bases.  The code reads an FCIDUMP
file, performs a GA simulation using diagonal-element information, and writes a
reordered FCIDUMP file using the best ordering found.  For details of the
algorithm, please see Section 2 of Reference [1].

## Quick Start

You can either directly execute `ga_driver.py` with input arguments, e.g.
```
# This is equivalent to `20site_Heisenberg_chain.py`.
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
```
or install the package with `pip` (see **Installation** below).

Once installed, you can import the package with `import GA_mod.run_GA as ga`.
Then, a GA simulation can be performed with `ga.perform_GA`.
Example GA scripts are located in the `examples/` directory.
During the course of a GA simulation, the latest population (set of orderings)
are stored in `current_pop.log`.
Below is a brief description of each:

- `fast_20site_Heisenberg_chain.py`<br/>
    This example runs a GA simulation for a 20-site NN Heisenberg chain
    (Subsection 4.1 of Reference [1]) using the simplified fitness function
    with the S-Ms mapping. Only ordering-dependent terms of the CSF energy
    are evaluated for fitness measure.

- `20site_Heisenberg_chain.py`<br/>
    This example is the same as `fast_20site_Heisenberg_chain`, but evaluates
    the full CSF energies, thus slower than `fast_20site_Heisenberg_chain`.

- `P-cluster_48i40.py`<br/>
    This example runs a GA simulation for the CAS(48e,40o) P-cluster consisting
    of 8 Fe sites (Subsection 5.2.1 of Reference [1]). The fitness function 
    evaluates the CSF energives over the 14 collinear VVS CSFs and maximizes the
    highest and the lowest CSF energies.

- `restart.py`<br/>
    This example runs a GA simulation starting from `NN_Heisenberg_20chain.pop`.

### Checkpoint feature

```bash
touch WRITE_CHECKPOINT
```

During a GA run, creating this file triggers writing an FCIDUMP file with the
current best ordering.

## Installation

```bash
cd GA4GUGA
pip install -e .
```

This installs the package in editable mode.

### Verify installation

```bash
python test_installation.py
```


## Fitness Functions

**GA4GUGA** supports fitness functions using CSF energies of various CSFs (see
Section 3 of Reference [1] for details) and using only S-Ms consistent CSF
energies (Subsection 3.1 of Reference [1]).

-  `MIN_MAX_DIFF`<br/>
    Maximizes the difference between the highest and the lowest CSF energies
    of a given set of CSFs.

-  `MAX_DIAG_ELEM`<br/>
    Maximizes the highest CSF energies of a given set of CSFs.

-  `MIN_DIAG_ELEM`<br/>
    Minimizes the lowest CSF energies of a given set of CSFs.

-  `DIAG_ELEM_SMS_MAPPING`<br/>
    Maximizes or minimizes (depending on the option) the energy of the S-Ms
    consistent CSF.

-  `FAST_DIAG_MIN_OSONLY`<br/>
    Minimizes the energy of the S-Ms consistent CSF, but only evaluates ordering
    -dependent terms in the CSF energy evaluation, thus faster than
    `DIAG_ELEM_SMS_MAPPING`.

## Input arguments

Arguments for the `run_GA.ga.perform_GA` funtion are listed below.
CSFs are assumed to be in the step-vector notation where $0$, $1$, $2$, and $3$
represent empty, increasing intermediate spin by 1/2, decreasing by 1/2,
and doubly occupied orbital couplings, respectively.

### Required arguments

- **fitness_function** `measure_fitness.FitnessFunction.FITNESS_FUNCTION_NAME`<br/>
    Chooses the fitness function for the GA simulation.
    Fitness functions are implemented in `GA_mod.measure_fitness.FitnessFunction`.

- **co_function** `crossover.CO_FUNCTION_NAME`<br/>
    Chooses the crossover operator for the GA simulation.
    Crossover operators are defined in `GA_mod.crossover`
    `order_co` is the order crossover operator.

- **num_chroms** $n$(int)<br/>
    Population size (number of chromosomes).

- **elite_size** $n$(int)<br/>
    Number of elite chromosomes retained each generation.

- **mutation_rates** $mlist$(list[float])<br/>
    A list of mutation rates. If multiple rates are provided, the algorithm can
    cycle through them when progress stalls. Even if a single mutation rate is
    provided, it has to be in list format (e.g., `[0.001]`).

- **restricted_ordering_len** $n$(int)<br/>
    Length of the effective (restricted) ordering that the GA optimizes before
    expansion (see **Restricting orderings** below for ordering expansion).

- **generations** $n$(int)<br/>
    Number of generations to run.

- **fcidump** `PATH_TO_FCIDUMP`<br/>
    Path to the source FCIDUMP file.

- **norb** $n$(int)<br/>
    Total number of orbitals. Must satisfy the consistency check:
    `norb = num_prefix + num_suffix + restricted_ordering_len * len(on_site_permutation)`

### Optional arguments

- **cluster_period** $n$(int)<br/>
    Period of cluster shuffling replacing the crossover step.
    The default value is 5.

- **stagnation_limit** $n$(int)<br/>
    If the best fitness score does not evolve for $n$ generations, the mutation
    rate is updated to the next value specified in `mutation_rates`.
    The default value is 100.

#### Fitness-function specific arguments

##### Diagonal elements of given CSFs (`MIN_MAX_DIFF`, `MAX_DIAG_ELEM`, and `MIN_DIAG_ELEM`)

- **csf_list** `csf_list`(list[list[int]])<br/>
    Specifies the list of CSFs that are used for fitness score evaluation.

##### S-Ms mapping fitness functions (`DIAG_ELEM_SMS_MAPPING` and `FAST_DIAG_MIN_OSONLY`)

- **sms_ref_csf** $csf$(list[int])<br/>
    Defines the S-Ms consistent CSF in the natural ordering.
    For example, `sms_ref_csf = [1, 2, 1, 2]` maps
    `{orbital:coupling} = {1:1, 2:2, 3:1, 4:2}`.

- **tMinimize** $l$(bool)<br/>
    Only required for `DIAG_ELEM_SMS_MAPPING` fitness function.
    This booliean determins whether the GA maximizes (if false) or minimizes
    (if true) the fitness score.

#### Restricting orderings
If `norb = 10`, `num_prefix = 3`, `num_suffix = 2`, `restricted_ordering_len = 5`, and `on_site_permutation = (1,)` (default), the GA is performed only for the middle set `{4,5,6,7,8}`. Orderings maintain the structure `(1,2,3, {4,5,6,7,8}, 9,10)`, e.g., `(1,2,3,4,6,7,5,8,9,10)` or `(1,2,3,6,8,7,4,5,9,10)`.
If `norb = 10`, `num_prefix = 2`, `num_suffix = 2`, `restricted_ordering_len = 3`, and `on_site_permutation = (2,1)`, then orderings maintain the structure `(1,2, {(4,3), (6,5), (8,7)}, 9,10)`, e.g., `(1,2,5,6,4,3,8,7,9,10)`.

- **on_site_permutation** $perm$(tuple[int])<br/>
    Pattern applied to each gene (number in the ordering) when expanding the
    restricted ordering to the full ordering. Default: `(1,)`.

- **num_prefix** $n$(int)<br/>
    Number of leading orbitals that remain fixed (not optimized).
    Default: `0`.

- **num_suffix** $n$(int)<br/>
    Number of trailing orbitals that remain fixed (not optimized).
    Default: `0`.

#### Restart

- **restart_filename** `pop_filename`<br/>
    If provided, restart the GA from a saved population.

#### Logging arguments

- **checkpoint_trigger** n(str)<br/>
    File name whose presence triggers writing an FCIDUMP file of the current
    best ordering during the run.
    Default: `"WRITE_CHECKPOINT"`.

## Citation

If you cite this package, please cite:

[1] M. Song and G. Li Manni, A Genetic Algorithm Approach for Compact Wave Function Representations in Spin-Adapted Bases, J. Chem. Theory Comput. **XXX**, XXXX (2025), DOI: [10.1021/acs.jctc.5c01264](https://doi.org/10.1021/acs.jctc.5c01264)
