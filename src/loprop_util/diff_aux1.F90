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

subroutine Diff_Aux1(nEPotPoints,EPCo,nB,OneFile)

use OneDat, only: sOpSiz
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: nEPotPoints
real(kind=wp), allocatable, intent(out) :: EPCo(:,:)
integer(kind=iwp), intent(in) :: nB
character(len=10), intent(in) :: OneFile
character(len=10) :: Label
integer(kind=iwp) :: i, iComp, iopt, irc, iSmLbl, Lu_One, maxCen, n_int(1)
real(kind=wp), allocatable :: idiot(:), Tmp(:,:)
integer(kind=iwp), external :: IsFreeUnit
#include "warnings.h"

! Open One-electron file.

irc = -1
Lu_One = 49
Lu_One = IsFreeUnit(Lu_One)
iopt = 0
call OpnOne(irc,iopt,OneFile,Lu_One)
if (irc /= 0) then
  write(u6,*)
  write(u6,*) 'ERROR! Could not open one-electron integral file.'
  call Quit(_RC_IO_ERROR_READ_)
end if

! Loop over all EF0, terminate when return-code is non-zero.

nEPotPoints = 0
maxCen = 99999
call mma_allocate(Tmp,3,maxCen,label='Temporary')
call mma_allocate(idiot,nB*(nB+1)/2+4,label='Idiot')
do i=1,maxCen
  write(Label,'(A3,I5)') 'EF0',i
  irc = -1
  iopt = ibset(0,sOpSiz)
  iSmLbl = 0
  iComp = 1
  call iRdOne(irc,iopt,label,iComp,n_Int,iSmLbl)
  if (irc /= 0) exit
  irc = -1
  iopt = 0
  iSmLbl = 0
  call RdOne(irc,iopt,label,iComp,idiot,iSmLbl)
  Tmp(:,i) = idiot(n_int(1)+1:n_int(1)+3)
  nEPotPoints = nEPotPoints+1
end do

! Put the coordinates and nuclear part in nice and tight vectors.

call mma_allocate(EPCo,3,nEPotPoints,label='PotPointCoord')
EPCo(:,:) = Tmp(:,1:nEPotPoints)

! Deallocate.

call mma_deallocate(Tmp)
call mma_deallocate(idiot)

return

end subroutine Diff_Aux1
