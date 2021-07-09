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

subroutine testq(nAt,nTs,VDer,Q,QTot)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs
real(kind=wp), intent(in) :: Q(2,*)
real(kind=wp), intent(_OUT_) :: VDer(nTs,*), QTot(*)
integer(kind=iwp) :: iAt, iCoord, idx, iTs, Lu
real(kind=wp) :: rsum

Lu = 1
call Molcas_open(Lu,'DerPt.dat')
!open(1,file='DerPot.dat',status='old',form='formatted')
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    do iTs=1,nTs
      read(1,*) VDer(iTs,idx)
    end do
  end do
end do
close(1)
do iAt=1,nAt
  do iCoord=1,3
    idx = 3*(iAt-1)+iCoord
    rsum = Zero
    do iTs=1,nts
      QTot(iTs) = Q(1,iTs)+Q(2,its)
      rsum = rsum+QTot(iTs)*VDer(iTs,idx)
    end do
    write(u6,'("Charges times VDer",i4,f20.12)') idx,rsum
  end do
end do

return

end subroutine testq
