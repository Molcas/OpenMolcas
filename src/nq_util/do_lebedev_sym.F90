!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine Do_Lebedev_Sym(L_Eff,mPt,ipR)

use stdalloc, only: mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L_Eff
integer(kind=iwp), intent(out) :: mPt, ipR
#include "WrkSpc.fh"
integer(kind=iwp) :: i, j, mPt_
real(kind=wp), allocatable :: R(:,:)
real(kind=wp), parameter :: Thr = 1.0e-16_wp
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Do_Lebedev(L_Eff,mPt,R)
    import :: wp, iwp
    implicit none
    integer(kind=iwp), intent(in) :: L_Eff
    integer(kind=iwp), intent(out) :: mPt
    real(kind=wp), allocatable, intent(out) :: R(:,:)
  end subroutine Do_Lebedev
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
call Do_Lebedev(L_Eff,mPt_,R)
mPt = 0
outer: do i=1,mPt_
  do j=1,i-1
    if (all(abs(R(1:3,j)+R(1:3,i)) < Thr)) then
      R(4,i) = Zero
      cycle outer
    end if
  end do
  mPt = mPt+1
end do outer

call GetMem('AngRW','Allo','Real',ipR,4*mPt)

j = 1
do i=1,mPt_
  if (R(4,i) /= Zero) then
    Work(ipR+(j-1)*4:ipR+j*4-1) = R(:,i)
    j = j+1
  end if
end do
call mma_deallocate(R)

end subroutine Do_Lebedev_Sym
