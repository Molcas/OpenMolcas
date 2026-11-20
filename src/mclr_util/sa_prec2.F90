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

subroutine SA_PREC2(rdia,S,CI,ENE)

use CI_Info, only: CI_St
use input_mclr, only: nCSF, nRoots
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
real(kind=wp), intent(in) :: rdia(nCSF(CI_St(iState)%STSym)), CI(nCSF(CI_St(iState)%STSym),nCSF(CI_St(iState)%STSym)), ENE
real(kind=wp), intent(out) :: S(nroots,nroots)
integer(kind=iwp) :: i, j, k, n
real(kind=wp) :: dnum

n = nCSF(CI_St(iState)%STSym)
S(:,:) = Zero
do i=1,nroots
  do j=1,nroots
    do k=1,n
      ! |H0-E|: keeps the preconditioner positive definite for excited
      ! roots (must be the same as in DMinvCI_SA)
      dnum = max(abs(rdia(k)-Ene),1.0e-16_wp)
      S(i,j) = S(i,j)+CI(k,i)*CI(k,j)/dnum
    end do
  end do
end do
call MatInvert(S,nroots)

end subroutine SA_PREC2
