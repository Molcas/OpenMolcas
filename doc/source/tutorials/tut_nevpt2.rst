.. index::
   single: Program; NEVPT2
   single: NEVPT2

.. _TUT\:sec\:nevpt2:

:program:`NEVPT2` --- :math:`n`-Electron Valence State Second-Order Perturbation Theory
=======================================================================================

NEVPT2 is a second-order perturbation theory with a CAS (or a CAS-like) reference wavefunction originally developed by Angeli et al. :cite:`Angeli_JChemPhys_Introduction_2001,Angeli_ChemPhysLett_Nelectron_2001,Angeli_JChemPhys_nelectron_2002,Angeli_JChemPhys_quasidegenerate_2004` In contrast to CASPT2, it uses a Dyall Hamiltonian :cite:`Dyall_JChemPhys_choice_1995` as the zeroth-order Hamiltonian and is therefore inherently free of intruder states and parameters such as the IPEA shift. NEVPT2 exists in two formulations -- the strongly- (SC-) and the partially-contracted NEVPT2 (PC-NEVPT2), which differ in the basis of the first-order wavefunction expansion.

The implementation in the :program:`NEVPT2` program is based on the original NEVPT2 implementation by Angeli et al. :cite:`Angeli_JChemPhys_nelectron_2002,Angeli_JChemPhys_quasidegenerate_2004`, with the implementation of the QCMaquis DMRG reference wave function and Cholesky decomposition for the two-electron integrals :cite:`Freitag_JChemTheoryComput_Multireference_2017`. For excited states both single-state and multi-state calculations with the QD-NEVPT2 approach :cite:`Angeli_JChemPhys_quasidegenerate_2004` are supported.

.. _TUT\:sec\:nevpt2_run:

Running a NEVPT2 calculation
----------------------------

Prior to running a NEVPT2 calculation, one must obtain a reference wavefunction with the :program:`RASSCF` or :program:`DMRGSCF` program and perform an integral transformation with the :program:`MOTRA` program.

Currently, the implementation supports **only** QCMaquis DMRG reference wavefunctions (support for CASSCF reference wavefunctions will be added in the near future). It is nevertheless possible to run NEVPT2 with a CASSCF reference wavefunction by performing a DMRG-CI calculation with a sufficiently large :math:`m` value using the CASSCF converged orbitals. For example, an :math:`m` value of 2000 recovers the exact CASCI energy up to :math:`5\times{}10^{-8}` a.u. for active spaces of up to 14 orbitals.

Below we show an example workflow of a NEVPT2 calculation. The input below is a calculation of the lowest singlet state of methylene with an active space of 6 electrons in 6 orbitals:

::

  &GATEWAY
    coord
    3
    CH2 Triplet coordinates in Angstrom
    C      0.000000  0.000000  0.000000
    H      0.000000  0.000000  1.077500
    H      0.784304  0.000000 -0.738832
    basis=cc-pVTZ
    Group=Nosym
    RICD
    CDTH=1.0E-7
  &SEWARD
  &DMRGSCF
    ActiveSpaceOptimizer=QCMaquis
    DMRGSettings
      max_bond_dimension=128
      nsweeps=5
    EndDMRGSettings
    OOptimizationSettings
      Spin=3
      Inactive=1
      Ras2=6
      NActEl=6,0,0
      NEVPT2Prep
    EndOOptimizationSettings
  &MOTRA
    Frozen=0
    CTOnly
    Kpq
    HDF5
  &NEVPT2

First, one performs a DMRG-SCF calculation with the keyword :kword:`NEVPT2Prep`, which enables the evaluation of the four-particle reduced density matrices (4-RDMs) (and, in case of multiple states, also transition three-particle density matrices (t-3RDMs)) required by NEVPT2.

Second, one must perform an integral transformation with the :program:`MOTRA` module. If no Cholesky decomposition or RICD is used in the calculation, the only mandatory keyword is :kword:`HDF5`, which enables the write-out of the transformed integrals in the HDF5 format required by the :program:`NEVPT2` module. If Cholesky decomposition is used, one additionally needs to add the keys :kword:`CTOnly` and :kword:`Kpq`. Cholesky decomposition is strongly recommended, as the integral transformation without Cholesky is several times slower and not supported in parallel.

Note that running with the Cholesky decomposed integrals currently does not support symmetry, and the support for frozen orbitals in :program:`MOTRA` with Cholesky is untested, hence also the keyword :kword:`Frozen=0` is recommended.

Finally, one calls the NEVPT2 module with :kword:`\&NEVPT2`. It has no mandatory options, but options described in the Users Guide can be specified.

.. _TUT\:sec\:nevpt2_distrdm:

Distributed RDM evaluation
--------------------------

The computational cost of the RDM evaluation grows as :math:`N^8` with the number of active orbitals, therefore the RDM evaluation for active spaces larger than 11-12 orbitals becomes prohibitively expensive. Therefore :program:`NEVPT2` distribution provides an (experimental) python utility :file:`jobmanager.py` for distributed massively parallel 4-RDM calculations. With distributed 4-RDM calculations, active spaces of up to 22 orbitals can be employed in DMRG-NEVPT2 calculations without any approximation to the 4-RDM.

:file:`jobmanager.py` splits the evaluation of the 4-RDM :math:`G_{ijklmnop}` into four-index subblocks with indices :math:`i,j,k,l`. Due to permutational symmetry and the properties of the creation and annihilation operators, :math:`i \ge j \ge k \ge l` and no more than two indexes are equal (pairwise equality :math:`i=j` and :math:`k=l` is allowed). The script prepares input files and, if requested, submits a separate job for each subblock, and merges the subblocks into the full matrix once the jobs are finished. The script is expected to be run on a head node of a distrubuted computing system with a batch system: `LSF <https://www.ibm.com/support/knowledgecenter/en/SSETD4/product_welcome_platform_lsf.html>`_ has been tested, but any batch system which supports the `DRMAA <http://www.drmaa.org/>`_ library, such as Slurm or PBS, should work. If no support for DRMAA is found, the script still may be used to prepare the input files for each subblock, which then may be submitted manually. Note that the DMRG-SCF/NEVPT2 calculation need not be performed on the same system as the 4-RDM evaluation.

How to run NEVPT2 calculations with distributed 4-RDM evaluation
................................................................

Prerequisites:

- Python :math:`\ge` 2.7.9 (3.x is also supported)

- (optional) DRMAA library compatible with your batch submission system, (e.g. `LSF-DRMAA <https://github.com/IBMSpectrumComputing/lsf-drmaa>`_)

- `Python DRMAA <https://github.com/drmaa-python/drmaa-python>`_

- (optional) GNU Parallel

If your system administrator has not set up DRMAA and Python DRMAA, you might need to download and install these libraries yourself. After the installation, the environment variable :variable:`DRMAA_LIBRARY_PATH` must be set to the path to :file:`libdrmaa.so` and, if Python does not find the DRMAA Python binding, also :variable:`PYTHONPATH` to the path of the Python DRMAA library.

Workflow:

- Run DMRGSCF and MOTRA calculations as shown above, but **omit** calling the :program:`NEVPT2` program. The :kword:`NEVPT2Prep` keyword in the :program:`DMRGSCF` section creates QCMaquis input templates and the MPS checkpoint files required for a later 4-RDM and/or t-3RDM evaluation.

- Copy the :file:`$MOLCAS/Tools/distributed-4rdm/prepare_rdm_template.sh` script to the |openmolcas| scratch directory and run it. The script will create subdirectories named :file:`4rdm-scratch.<state>` for each state. If you wish to perform the 4-RDM evaluation on a different machine (e.g. a cluster), copy the subdirectory for each state to that machine. If you do not wish to evaluate the 4-RDM for all states, pass the list of desired states as parameters to the :file:`prepare_rdm_template.sh` script. For example, :file:`./prepare_rdm_template.sh 0 1 2` will create the scratch directories for states from 0 to 2 (note that QCMaquis starts counting states with 0).

- **If you have installed and working DRMAA setup:** For each state, change to the :file:`4rdm-scratch.<state>` subdirectory and run ::

     nohup jobmanager.py &

  (Login to the machine where you evaluate the 4-RDM before if you wish to run the evaluation on a different machine.) This will create a subdirectory for each batch job (corresponding for each four-index 4-RDM subblock) and submit the jobs. The script will stay in the background until all the jobs have completed.
  The script also accepts the following job-specific options:

  - :command:`-t HH:MM:SS`: set the maximum walltime per job. Default is 24h.
  - :command:`-n NCPU`: run each job in an SMP parallelised fashion and set the number of CPU cores per job. Default is 1 core. For large active spaces, it is recommended to use several cores (e.g. 16 or 24, or as much as is available on a single node on your cluster).

- If you **do not** have DRMAA installed and working, run the :file:`jobmanager.py` script with the :command:`-n` option: ::

    jobmanager.py -n

  This will create subfolders for each 4-RDM block and prepare all the necessary input scripts, but will not submit them to the batch system. Now you may manually submit the scripts from the subfolders :file:`parts/part-*`.

- If you ran the distributed 4-RDM calculation on a different machine, copy the :file:`4rdm-scratch.<state>` back to |openmolcas| :file:`$WorkDir`.

- Create an input file with the input to the :program:`NEVPT2` program and run it. The keyword :kword:`DISTributedRDM` followed by the path to :file:`4rdm-scratch.<state>` folders (in our case, :file:`$WorkDir`) is **mandatory**.

Troubleshooting
...............

The :file:`jobmanager.py` script is experimental, and also batch jobs in queuing systems are prone to crash, therefore we provide a mechanism to identify and restart the crashed batch jobs. The NEVPT2 program will check if the 4-RDM calculation has been finished correctly. If some 4-RDM values are missing, the NEVPT2 program will stop with an error. In this case several options are available:

- **If DRMAA has been used:** if the :file:`jobmanager.py` finishes without errors, it will produce two files, :file:`successlist` and :file:`faillist` with the list of successful and failed batch jobs, respectively. In this case, the failed jobs may be restarted using the restart mode of :file:`jobmanager.py`, which is invoked with ::

    nohup jobmanager.py -r successlist faillist &

  If the :file:`jobmanager.py` finishes with an error, the :file:`successlist` and :file:`faillist` will be either nonexistent or empty. Note that this does NOT necessarily mean that the jobs have failed: in our tests, certain configurations of the queuing system may lead to the crash of the :file:`jobmanager.py` script after the successful completion of the jobs.

- **If DRMAA has not been used and the script was run with the -n switch**: in this case the user is advised to check manually the subfolders for each 4-RDM subblock for the existence of :file:`$Project.results_state.X.h5` files. The files should exist and the command

  .. compound::

     ::

       h5dump $Project.results_$state.X.h5 | grep fourpt

     should not yield an empty result -- otherwise the corresponding calculation should be rerun.

- Finally, if :program:`NEVPT2` is started with the :kword:`DISTributedRDM` keyword, it will check the number of evaluated 4-RDM elements. If the number of evaluated elements is different from its expected value, the program will exit with an error.

Transition 3-RDM distributed calculations
.........................................

:file:`jobmanager.py` also supports distributed calculations of t-3RDMs (required for multi-state QD-NEVPT2). The split evaluation is similar to that of the 4-RDMs, and the workflow above can be followed with the following differences:

- The t-3RDM evaluation requires two states instead of one. Run the :file:`prepare_rdm_template.sh` script with the :command:`-3` parameter.

- Launch the :file:`jobmanager.py` script with the :command:`-3` parameter.
