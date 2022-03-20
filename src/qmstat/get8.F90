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
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Ract, Etot
#include "maxi.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, iCT, iDisk, j
real(kind=wp) :: Esub, Gamold, GaOld
character(len=200) :: Head

iDisk = 0
call DaName(iLuStIn,StFilIn)
call WrRdSim(iLuStIn,2,iDisk,iTcSim,64,Etot,Ract,nPart,Gamold,GaOld,Esub)
iDisk = iTcSim(1)

! In this loop we read the coordinates. The construction of Cordst
! makes this loop necessary. Maybe we should consider going to
! dynamic allocation.

do i=1,3
  call GetMem('CTemp','Allo','Real',iCT,nPart*nCent)
  call dDaFile(iLuStIn,2,Work(iCT),nPart*nCent,iDisk)
  do j=1,nCent*nPart
    Cordst(j,i) = Work(iCT+j-1)
  end do
  call GetMem('CTemp','Free','Real',iCT,nPart*nCent)
  iDisk = iTcSim(i+1)
end do
call DaClos(iLuStIn)

! If requested, print initial coordinates.

if (iPrint >= 10) then
  write(Head,*) 'Coordinates read from startfile.'
  call Cooout(Head,Cordst,nPart,nCent)
end if

return

end subroutine Get8
