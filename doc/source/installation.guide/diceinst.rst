.. _sec\:dice_installation:

Installation of Dice--|molcas| interface for SHCI calculations
=================================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none


The Dice--|molcas| interface allows one to use semistochastic heat-bath configuration interaction (SHCI)
implemented in Dice as an FCI solver in CASSCF calculations, referred to as SHCI-CASSCF :cite:`Phung2016,Wouters2016`.
A large active space, up to around 100 active orbitals, can be calculated with SHCI-CASSCF.
Currently, the interface supports ground state SHCI-CASSCF calculations.

The interface requires the Dice 1.0 binary (https://github.com/sanshar/Dice).
For installation of Dice, consult https://sanshar.github.io/Dice/installation.html.
The interface supports both parallel Dice and Molcas.

The Dice--|molcas| interface is built by activating:

::

  -D DICE=ON

Before runing SHCI-CASSCF calculations with the Dice-|molcas| interface, make sure to increase stack size;
and export the Dice binary and all the required libraries for Dice.

::

  ulimit -s unlimited
  export PATH=/path/to/dice/binary:$PATH

To run parallel Dice, export the environment variable MOLCAS_DICE, for example when running on 16 nodes use:

::

  export MOLCAS_DICE=16

Verify the installation:

::

  molcas verify extra:852-854
