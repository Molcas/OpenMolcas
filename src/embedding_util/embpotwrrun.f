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
      subroutine embPotWrRun

      implicit none

      ! Holds the data which is written in this subroutine
#include "embpotdata.fh"

      integer ioUnit, isFreeUnit

      ! Open the file
      ioUnit = isFreeUnit(1)
      call molcas_open(ioUnit, 'EMBPOT_RUNFILE')

      write(ioUnit,'(L1)') embPotInBasis
      write(ioUnit,'(A256)') embPotPath
      write(ioUnit,'(L1)') outGridPathGiven
      write(ioUnit,'(A256)') outGridPath
      write(ioUnit,'(L1)') embWriteDens
      write(ioUnit,'(A256)') embOutDensPath
      write(ioUnit,'(L1)') embWriteEsp
      write(ioUnit,'(A256)') embOutEspPath
      write(ioUnit,'(L1)') embWriteGrad
      write(ioUnit,'(A256)') embOutGradPath
      write(ioUnit,'(L1)') embWriteHess
      write(ioUnit,'(A256)') embOutHessPath

      close(ioUnit)

      return

      end
