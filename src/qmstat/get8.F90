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

! GET NUMBERS FROM STARTFILE.
subroutine Get8(Ract,Etot)

use qmstat_global, only: Cordst, iLuStIn, iPrint, iTcSim, nCent, nPart, StFilIn
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(out) :: Ract, Etot
integer(kind=iwp) :: i, iDisk
real(kind=wp) :: Esub, Gamold, GaOld
character(len=200) :: Head
real(kind=wp), allocatable :: CT(:)

iDisk = 0
call DaName(iLuStIn,StFilIn)
call WrRdSim(iLuStIn,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,GaOld,Esub)
iDisk = iTcSim(1)

! In this loop we read the coordinates. The construction of Cordst
! makes this loop necessary. Maybe we should consider going to
! dynamic allocation.

call mma_allocate(CT,nPart*nCent,label='CTemp')
do i=1,3
  call dDaFile(iLuStIn,2,CT,nPart*nCent,iDisk)
  Cordst(i,:) = CT
  iDisk = iTcSim(i+1)
end do
call mma_deallocate(CT)
call DaClos(iLuStIn)

! If requested, print initial coordinates.

if (iPrint >= 10) then
  write(Head,*) 'Coordinates read from startfile.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine Get8
