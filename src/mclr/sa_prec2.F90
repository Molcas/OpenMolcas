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

use input_mclr, only: nCSF, nRoots, State_Sym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: rdia(nCSF(State_Sym)), CI(nCSF(State_Sym),nCSF(State_Sym)), ENE
real(kind=wp), intent(out) :: S(nroots,nroots)
integer(kind=iwp) :: i, j, k, n
real(kind=wp) :: dnum

n = nCSF(State_Sym)
S(:,:) = Zero
do i=1,nroots
  do j=1,nroots
    do k=1,n
      dnum = rdia(k)-Ene
      dnum = sign(max(abs(dnum),1.0e-16_wp),dnum)
      S(i,j) = S(i,j)+CI(k,i)*CI(k,j)/dnum
    end do
  end do
end do
call MatInvert(S,nroots)

end subroutine SA_PREC2
