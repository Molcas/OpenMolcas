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

subroutine embPotWrRun()

use Embedding_Global, only: embOutDensPath, embOutEspPath, embOutGradPath, embOutHessPath, embPotInBasis, embPotPath, &
                            embWriteDens, embWriteEsp, embWriteGrad, embWriteHess, outGridPath, outGridPathGiven
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ioUnit
integer(kind=iwp), external :: isFreeUnit

! Open the file
ioUnit = isFreeUnit(1)
call molcas_open(ioUnit,'EMBPOT_RUNFILE')

write(ioUnit,'(l1)') embPotInBasis
write(ioUnit,'(a256)') embPotPath
write(ioUnit,'(l1)') outGridPathGiven
write(ioUnit,'(a256)') outGridPath
write(ioUnit,'(l1)') embWriteDens
write(ioUnit,'(a256)') embOutDensPath
write(ioUnit,'(l1)') embWriteEsp
write(ioUnit,'(a256)') embOutEspPath
write(ioUnit,'(l1)') embWriteGrad
write(ioUnit,'(a256)') embOutGradPath
write(ioUnit,'(l1)') embWriteHess
write(ioUnit,'(a256)') embOutHessPath

close(ioUnit)

return

end subroutine embPotWrRun
