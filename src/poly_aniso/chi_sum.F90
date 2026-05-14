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
!    Xex -- susceptibility tensor arising form the exchange states only, Real(kind=wp) :: , (3,3) array, input
!    Zex -- statistical sum according to Boltzmann distribution law of the exchange states, Real(kind=wp) :: , input
!    XL  -- susceptibility tensor arising from LOCAL states (all), Real(kind=wp) :: , (N,3,3) array, input
!    ZL  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (all),
!           Real(kind=wp) :: , (N) array, input
!    XR  -- susceptibility tensor arising from LOCAL states (exchange only), Real(kind=wp) :: , (N,3,3) array, input
!    ZR  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (exchange only),
!           Real(kind=wp) :: , (N) array, input
!   iopt -- option allowing to choose the desired formula, (Integer, input):
!           iopt=1  =>  formula for weak exchange limit ( new derivation)
!           iopt=2  =>  formula for strong exchange limit
!           iopt=3  =>  formula for strong exchange limit
!    X   -- total susceptibility, Real(kind=wp) :: , (3,3) array, output
!    Z   -- total statistical sum according to Boltzmann distribution, Real(kind=wp) :: , output
!---------

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, iopt
real(kind=wp), intent(in) :: Xex(3,3), Zex, XL(N,3,3), ZL(N), XR(N,3,3), ZR(N)
real(kind=wp), intent(out) :: X(3,3), Z
integer(kind=iwp) :: i
real(kind=wp) :: XLT(3,3), XRT(3,3), ZLT, ZRT

XRT(:,:) = Zero
XLT(:,:) = Zero
! compute the total ZT
if (iopt == 1) then
  ! my formula (simple):
  ZLT = product(ZL)
  ZRT = product(ZR)
  Z = Zex+ZLT-ZRT
  do i=1,N
    XLT(:,:) = XLT(:,:)+XL(i,:,:)
    XRT(:,:) = XRT(:,:)+XR(i,:,:)
  end do
  X(:,:) = Xex(:,:)+XLT(:,:)-XRT(:,:)

else if (iopt == 2) then
  ! thesis formula:
  ZLT = product(ZL)
  ZRT = product(ZR)
  Z = Zex+ZLT-ZRT
  do i=1,N
    XLT(:,:) = XLT(:,:)+XL(i,:,:)*ZLT
    XRT(:,:) = XRT(:,:)+XR(i,:,:)*ZRT
  end do
  X(:,:) = (Xex(:,:)*Zex+XLT(:,:)-XRT(:,:))/Z

else if (iopt == 3) then
  ! weird formula as implemented in some version of the code, e.g. in the 2013 version:
  ZLT = product(ZL)
  ZRT = product(ZR)
  ZLT = Zex*ZLT/ZRT
  ZRT = Zex
  Z = ZLT-ZRT+Zex
  do i=1,N
    XLT(:,:) = XLT(:,:)+XL(i,:,:)
    XRT(:,:) = XRT(:,:)+XR(i,:,:)
  end do
  X(:,:) = (XLT(:,:)*ZLT-XRT(:,:)*ZRT+Xex(:,:)*Zex)/Z

else

  Z = Zero
  X(:,:) = Zero
  write(u6,'(A)') 'chi_sum: IOPT parameter out of range'
  write(u6,'(A,i8)') 'IOPT = ',IOPT
  return

end if

return

end subroutine chi_sum
