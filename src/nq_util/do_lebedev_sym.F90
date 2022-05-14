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

! This subroutine should be in a module, to avoid explicit interfaces
#ifdef _IN_MODULE_

subroutine Do_Lebedev_Sym(L_Eff,mPt,R)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L_Eff
integer(kind=iwp), intent(out) :: mPt
real(kind=wp), allocatable, intent(out) :: R(:,:)
integer(kind=iwp) :: i, j, mPt_f
real(kind=wp), allocatable :: R_f(:,:)
real(kind=wp), parameter :: Thr = 1.0e-16_wp

!                                                                      *
!***********************************************************************
!                                                                      *
call Do_Lebedev(L_Eff,mPt_f,R_f)
mPt = 0
outer: do i=1,mPt_f
  do j=1,i-1
    if (all(abs(R_f(1:3,j)+R_f(1:3,i)) < Thr)) then
      R_f(4,i) = Zero
      cycle outer
    end if
  end do
  mPt = mPt+1
end do outer

call mma_allocate(R,4,mPt,label='R')

j = 1
do i=1,mPt_f
  if (R_f(4,i) /= Zero) then
    R(:,j) = R_f(:,i)
    j = j+1
  end if
end do
call mma_deallocate(R_f)

end subroutine Do_Lebedev_Sym

#endif
