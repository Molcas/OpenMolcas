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

use input_mclr, only: nRoots, nCSF, State_Sym
use Constants, only: Zero
use Definitions, only: wp

implicit none
real*8 rdia(*), S(nroots,nroots), CI(*)
real*8 ENE
integer i, j, k
real*8 dnum

do i=0,nroots-1
  do j=0,nroots-1
    S(i+1,j+1) = Zero
    do k=1,ncsf(State_Sym)
      dnum = rdia(k)-Ene
      dnum = sign(max(abs(dnum),1.0e-16_wp),dnum)
      S(i+1,j+1) = S(i+1,j+1)+CI(i*ncsf(State_Sym)+k)*CI(j*ncsf(State_Sym)+k)/dnum
    end do
  end do
end do
call MatInvert(S,nroots)

end subroutine SA_PREC2
