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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine HighFm(F,T,m,n)
!***********************************************************************
!  Object: to compute the auxiliary function for orders which we do    *
!          not use Shavitt's method of tables.                         *
!                                                                      *
! Called from: Auxil                                                   *
!                                                                      *
! Calling    : None                                                    *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!***********************************************************************

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: m, n
real(kind=wp), intent(out) :: F(n)
real(kind=wp), intent(in) :: T(n)
integer(kind=iwp) :: i, ii, k
real(kind=wp) :: FValue, Gamma2, gTmp, Sum0, Sum1, Term, TMax, TNew

! Find T for which the asympotic formula can be used

Tmax = 50.0_wp
do
  gTmp = Gamma2(m,Tmax)
  i = 1
  ii = 2*m-1
  sum1 = One
  sum0 = One
  do
    Sum1 = Sum1*real(ii,kind=wp)/(Two*Tmax)
    Sum0 = Sum0+Sum1
    ii = ii-2
    i = i+1
    if ((i < m) .or. (Sum1/sum0 <= 1.0e-11_wp)) exit
  end do
  Tnew = log(sum0/(2.0e-16_wp*Tmax*gTmp))
  if (abs(Tnew-Tmax) < 1.0e-9_wp) exit
  Tmax = Tnew
end do

! Compute the auxiliary functions

Tmax = Tnew
do k=1,n
  if (T(k) < Tmax) then
    Fvalue = Zero
    i = 0
    Term = One
    do
      Term = Term/real(2*(m+i)+1,kind=wp)
      Fvalue = Fvalue+Term
      i = i+1
      Term = Term*(Two*T(k))
      if (abs(Term/Fvalue) <= 1.0e-18_wp) exit
    end do
    F(k) = exp(-T(k))*Fvalue
  else
    F(k) = Gamma2(m,T(k))
  end if
end do

return

end subroutine HighFm
