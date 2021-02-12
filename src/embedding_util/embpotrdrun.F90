!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine embPotRdRun()

use Embedding_Global, only: embOutDensPath, embOutEspPath, embOutGradPath, embOutHessPath, embPotInBasis, embPotPath, &
                            embWriteDens, embWriteEsp, embWriteGrad, embWriteHess, outGridPath, outGridPathGiven
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ioUnit
integer(kind=iwp), external :: isFreeUnit

! Open the file
ioUnit = isFreeUnit(1)
call molcas_open(ioUnit,'EMBPOT_RUNFILE')

read(ioUnit,'(l1)') embPotInBasis
read(ioUnit,'(a256)') embPotPath
read(ioUnit,'(l1)') outGridPathGiven
read(ioUnit,'(a256)') outGridPath
read(ioUnit,'(l1)') embWriteDens
read(ioUnit,'(a256)') embOutDensPath
read(ioUnit,'(l1)') embWriteEsp
read(ioUnit,'(a256)') embOutEspPath
read(ioUnit,'(l1)') embWriteGrad
read(ioUnit,'(a256)') embOutGradPath
read(ioUnit,'(l1)') embWriteHess
read(ioUnit,'(a256)') embOutHessPath

close(ioUnit)

return

end subroutine embPotRdRun
