************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_SetDecAlg_Def()
C
C     Set default decomposition algorithm.
C
      Implicit None
#include "cholesky.fh"

CSVC This is what CHO_DECALG means according to TBP:
C     1 - serial one-step (the original algorithm)
C     2 - serial two-step (the original two-step algorithm)
C     3 - naive (mainly for test purposes; should never be used)
C     4 - parallel one-step (parallel modifications of the original one-step, can also run in serial - but gives a different result than the original serial one-step)
C     5 - parallel two-step (a two-step algorithm different from the serial one; can also run in serial)
C     6 - parallel naive (again, mainly for testing)

CSVC: the default used to be 5 (parallel two-step), but because of a problem
C with DGA where it hangs when using a largish number of processes, the default
C has been changed to 4, i.e. the parallel one-step algorithm. The default
C should be changed back to 5 if and when DGA and/or Cholesky have been adapted.

      Cho_DecAlg_Def=4

      End
