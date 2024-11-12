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

subroutine Msum(N,Mex,Zex,ML,ZL,MR,ZR,iopt,M,Z)
! computes the total M, provided all input values are provided
! according to the desired partition scheme  (iopt = 1 or 2 )
! all M vectors must be given in the general coordinate system
!
! definition of the variables:
!     N  -- total number of magnetic sites, Integer, input
!    Mex -- magnetisation vector arising from the exchange states only, Real(kind=wp) ::, (3) array, input
!    Zex -- statistical sum according to Boltzmann distribution law of the exchange states, Real(kind=wp) ::, input
!    ML  -- magnetisation vector arising from LOCAL states (all), Real(kind=wp) ::, (N,3) array, input
!    ZL  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (all),
!           Real(kind=wp) ::, (N) array, input
!    MR  -- magnetisation vector arising from LOCAL states (exchange only), Real(kind=wp) ::, (N,3) array, input
!    ZR  -- statistical sum according to Boltzmann distribution law arising from LOCAL states (exchange only),
!           Real(kind=wp) ::, (N) array, input
!   iopt -- option allowing to choose the desired formula, (Integer, input):
!           iopt=1  =>  formula for weak exchange limit ( new derivation)
!           iopt=2  =>  formula for strong exchange limit ( simple sumation of moments),size consistent;
!           iopt=3  =>  formula for strong exchange limit ( not to be used...)
!    M   -- total magnetisation, Real(kind=wp) ::, (3) array, output
!    Z   -- total statistical sum according to Boltzmann distribution, Real(kind=wp) ::, output
!---------

use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, iopt
real(kind=wp), intent(in) :: Mex(3), Zex, ML(N,3), ZL(N), MR(N,3), ZR(N)
real(kind=wp), intent(out) :: M(3), Z
integer(kind=iwp) :: i
real(kind=wp) :: MLT(3), MRT(3), ZLT, ZRT

MLT(:) = Zero
MRT(:) = Zero

if (iopt == 1) then
  ! compute the total ZT
  ! my formula (simple):
  ZLT = product(ZL)
  ZRT = product(ZR)
  Z = Zex+ZLT-ZRT
  do i=1,N
    MLT(:) = MLT(:)+ML(i,:)
    MRT(:) = MRT(:)+MR(i,:)
  end do
  M(:) = Mex(:)+MLT(:)-MRT(:)
# ifdef _DEBUGPRINT_
  write(u6,'(A,3F10.6)') 'Contribution from exchange states',(Mex(i),i=1,3)
  write(u6,'(A,3F10.6)') 'Contribution from excited states ',(MLT(i)-MRT(i),i=1,3)
# endif

else if (iopt == 2) then
  ! thesis formula:
  ZLT = product(ZL)
  ZRT = product(ZR)
  Z = Zex+ZLT-ZRT
  do i=1,N
    MLT(:) = MLT(:)+ML(i,:)*ZLT
    MRT(:) = MRT(:)+MR(i,:)*ZRT
  end do
  M(:) = (Mex(:)*Zex+MLT(:)-MRT(:))/Z

else

  Z = Zero
  M(:) = Zero
  write(u6,'(A)') 'chi_sum: IOPT parameter out of range'
  write(u6,'(A,i8)') 'IOPT = ',IOPT

end if

return

end subroutine Msum
