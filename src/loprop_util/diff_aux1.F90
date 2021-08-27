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

subroutine Diff_Aux1(nEPotPoints,ipEPCo,nB,OneFile)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nEPotPoints, ipEPCo
integer(kind=iwp), intent(in) :: nB
character(len=10), intent(in) :: OneFile
character(len=10) :: Label
integer(kind=iwp) :: i, idiot, iopt, irc, iSmLbl, iTmp, Lu_One, maxCen, n_int(1)
integer(kind=iwp), external :: IsFreeUnit
#include "WrkSpc.fh"
#include "warnings.h"

! Open One-electron file.

irc = -1
Lu_One = 49
Lu_One = IsFreeUnit(Lu_One)
call OpnOne(irc,0,OneFile,Lu_One)
if (irc /= 0) then
  write(u6,*)
  write(u6,*) 'ERROR! Could not open one-electron integral file.'
  call Quit(_RC_IO_ERROR_READ_)
end if

! Loop over all EF0, terminate when return-code is non-zero.

nEPotPoints = 0
maxCen = 99999
call GetMem('Temporary','Allo','Real',iTmp,maxCen*3)
call GetMem('Idiot','Allo','Real',idiot,nB*(nB+1)/2+4)
do i=1,maxCen
  write(Label,'(A3,I5)') 'EF0',i
  irc = -1
  iopt = 1
  iSmLbl = 0
  call iRdOne(irc,iopt,label,1,n_Int,iSmLbl)
  if (irc /= 0) exit
  irc = -1
  iopt = 0
  iSmLbl = 0
  call RdOne(irc,iopt,label,1,Work(idiot),iSmLbl)
  Work(iTmp+(i-1)*3+0) = Work(idiot+n_int(1)+0)
  Work(iTmp+(i-1)*3+1) = Work(idiot+n_int(1)+1)
  Work(iTmp+(i-1)*3+2) = Work(idiot+n_int(1)+2)
  nEPotPoints = nEPotPoints+1
end do

! Put the coordinates and nuclear part in nice and tight vectors.

call GetMem('PotPointCoord','Allo','Real',ipEPCo,3*nEPotPoints)
call dcopy_(3*nEPotPoints,Work(iTmp),1,Work(ipEPCo),1)

! Deallocate.

call GetMem('Temporary','Free','Real',iTmp,maxCen*3)
call GetMem('Idiot','Free','Real',idiot,nB*(nB+1)/2+4)

return

end subroutine Diff_Aux1
