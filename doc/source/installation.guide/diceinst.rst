.. _sec\:dice_installation:

Installation of Dice--|molcas| interface for HCI calculations
=============================================================

.. only:: html

  .. contents::
     :local:
     :backlinks: none


The Dice--|molcas| interface allows one to use heat-bath configuration interaction (HCI)
implemented in Dice as an FCI solver in CASSCF calculations, referred to as HCI-CASSCF :cite:`Sharma2017,Holmes2016`.
A large active space, up to around 100 active orbitals, can be calculated with HCI-CASSCF.
Currently, the interface supports ground state HCI-CASSCF calculations.

The interface requires the Dice 1.0 binary (https://github.com/sanshar/Dice).
For installation of Dice, consult https://sanshar.github.io/Dice/installation.html.
The interface supports both parallel Dice and |molcas|.

The Dice--|molcas| interface is built by activating in CMake:

::

  -D DICE=ON

Before runing HCI-CASSCF calculations with the Dice--|molcas| interface, make sure to increase stack size;
and export the Dice binary and all the required libraries for Dice.

::

  ulimit -s unlimited
  export PATH=/path/to/dice/binary:$PATH

To run parallel Dice, export the environment variable :variable:`MOLCAS_DICE`, for example when running on 16 nodes use:

::

  export MOLCAS_DICE=16

Verify the installation:

::

  molcas verify .all -w dice
