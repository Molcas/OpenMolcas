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

subroutine chi_sum(N,Xex,Zex,XL,ZL,XR,ZR,iopt,X,Z)
! computes the total CHI, provided all input values are provided
! according to the desired partition scheme  (iopt = 1 or 2 )
! all X tensors must be given in the general coordinate system
!
! definition of the variables:
!     N  -- total number of magnetic sites, Integer, input
!    Xex -- susceptibility tensor arising form the exchange states only, Real(kind=8) :: , (3,3) array, input
!    Zex -- statistical sum according to Boltzmann distribution law of the exchange states, Real(kind=8) :: , input
!    XL  -- susceptibility tensor arising from LOCAL states (all), Real(kind=8) :: , (N,3,3) array, input
!    ZL  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (all),
!           Real(kind=8) :: , (N) array, input
!    XR  -- susceptibility tensor arising from LOCAL states (exchange only), Real(kind=8) :: , (N,3,3) array, input
!    ZR  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (exchange only),
!           Real(kind=8) :: , (N) array, input
!   iopt -- option allowing to choose the desired formula, (Integer, input):
!           iopt=1  =>  formula for weak exchange limit ( new derivation)
!           iopt=2  =>  formula for strong exchange limit
!           iopt=3  =>  formula for strong exchange limit
!    X   -- total susceptibility, Real(kind=8) :: , (3,3) array, output
!    Z   -- total statistical sum according to Boltzmann distribution, Real(kind=8) :: , output
!---------

use Constants, only: Zero, One
use Definitions, only: u6

implicit none
integer, intent(in) :: N, iopt
real(kind=8), intent(in) :: Xex(3,3), Zex
real(kind=8), intent(in) :: XL(N,3,3), ZL(N)
real(kind=8), intent(in) :: XR(N,3,3), ZR(N)
real(kind=8), intent(out) :: X(3,3), Z
! local variables
integer :: i, ic, jc
real(kind=8) :: ZLT, ZRT, XLT(3,3), XRT(3,3)

ZLT = One
ZRT = One
Z = Zero
call dcopy_(3*3,[Zero],0,XRT,1)
call dcopy_(3*3,[Zero],0,XLT,1)
call dcopy_(3*3,[Zero],0,X,1)
! compute the total ZT
if (iopt == 1) then
  ! my formula (simple):
  do i=1,N
    ZLT = ZLT*ZL(i)
    ZRT = ZRT*ZR(i)
  end do
  Z = Zex+ZLT-ZRT
  do ic=1,3
    do jc=1,3
      do i=1,N
        XLT(ic,jc) = XLT(ic,jc)+XL(i,ic,jc)
        XRT(ic,jc) = XRT(ic,jc)+XR(i,ic,jc)
      end do
      X(ic,jc) = Xex(ic,jc)+XLT(ic,jc)-XRT(ic,jc)
    end do
  end do

else if (iopt == 2) then
  ! "thesis formula:"
  do i=1,N
    ZLT = ZLT*ZL(i)
    ZRT = ZRT*ZR(i)
  end do
  Z = Zex+ZLT-ZRT
  do ic=1,3
    do jc=1,3
      do i=1,N
        XLT(ic,jc) = XLT(ic,jc)+XL(i,ic,jc)*ZLT
        XRT(ic,jc) = XRT(ic,jc)+XR(i,ic,jc)*ZRT
      end do
      X(ic,jc) = (Xex(ic,jc)*Zex+XLT(ic,jc)-XRT(ic,jc))/Z
    end do
  end do

else if (iopt == 3) then
  ! "weird formula as implemented in some version of the
  ! code, e.g. in the 2013 version:"
  do i=1,N
    ZLT = ZLT*ZL(i)
    ZRT = ZRT*ZR(i)
  end do
  ZLT = Zex*ZLT/ZRT
  ZRT = Zex
  Z = ZLT-ZRT+Zex
  do ic=1,3
    do jc=1,3
      do i=1,N
        XLT(ic,jc) = XLT(ic,jc)+XL(i,ic,jc)
        XRT(ic,jc) = XRT(ic,jc)+XR(i,ic,jc)
      end do
      X(ic,jc) = (XLT(ic,jc)*ZLT-XRT(ic,jc)*ZRT+Xex(ic,jc)*Zex)/Z
    end do
  end do

else

  write(u6,'(A)') 'chi_sum: IOPT parameter out of range'
  write(u6,'(A,i8)') 'IOPT = ',IOPT
  return

end if

return

end subroutine chi_sum
