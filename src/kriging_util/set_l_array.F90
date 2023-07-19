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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine set_l_Array(Array_l,nInter,BaseLine,Hessian,HDiag)

use Constants, only: Three, Five
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(out) :: Array_l(nInter)
real(kind=wp), intent(in) :: BaseLine
real(kind=wp), intent(inout), optional :: Hessian(nInter,nInter), HDiag(nInter)
integer(kind=iwp) :: i
real(kind=wp) :: Hss
real(kind=wp), parameter :: H_min = 0.025_wp

!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
!call RecPrt('set_l_Array: Hessian',' ',Hessian,nInter,nInter)
!write(u6,*) 'BaseLine=',BaseLine
#endif

! Gives a Kriging Hessian for a single point of Kriging with
! a diagonal which is identical to the diagonal values of
! the HMF ad hoc Hessian.

if (present(Hessian)) then
  do i=1,nInter

    ! Make sure that the characteristic length is not too long.

    Hss = max(abs(Hessian(i,i)),H_min)
    Hessian(i,i) = Hss
    Array_l(i) = sqrt(Five/Three*BaseLine/Hss)

  end do
else
  do i=1,nInter

    ! Make sure that the characteristic length is not too long.

    Hss = max(abs(HDiag(i)),H_min)
    HDiag(i) = Hss
    Array_l(i) = sqrt(Five/Three*BaseLine/Hss)

  end do
end if

#ifdef _DEBUGPRINT_
!write(u6,*) 'H_min=',H_Min
!call RecPrt('Array_l',' ',Array_l,1,nInter)
#endif

return

end subroutine set_l_Array
