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
      subroutine embPotRdRun

      implicit none

      ! Holds the data which is read in in this subroutine
#include "embpotdata.fh"

      integer ioUnit, isFreeUnit

      ! Open the file
      ioUnit = isFreeUnit(1)
      call molcas_open(ioUnit, 'EMBPOT_RUNFILE')

      read(ioUnit,'(L1)') embPotInBasis
      read(ioUnit,'(A256)') embPotPath
      read(ioUnit,'(L1)') outGridPathGiven
      read(ioUnit,'(A256)') outGridPath
      read(ioUnit,'(L1)') embWriteDens
      read(ioUnit,'(A256)') embOutDensPath
      read(ioUnit,'(L1)') embWriteEsp
      read(ioUnit,'(A256)') embOutEspPath
      read(ioUnit,'(L1)') embWriteGrad
      read(ioUnit,'(A256)') embOutGradPath
      read(ioUnit,'(L1)') embWriteHess
      read(ioUnit,'(A256)') embOutHessPath

      close(ioUnit)

      return

      end
