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

subroutine FULLTRNSF(NP,NW,NB,CMOBLK,NJ,BUF_HT,BUF_FT)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NP, NW, NB, NJ
real(kind=wp), intent(in) :: CMOBLK(NB,NP), BUF_HT(NW*NJ,NB)
real(kind=wp), intent(out) :: BUF_FT(NP,NW*NJ)

! In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)

! Compute fully transformed Cholesky vector buffer:
!  FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
!do J=1,NJ
!  call DGEMM_('T','T',NP,NW,NB,One,CMOBLK,NB,BUF_HT(1,J,1),NW*NJ,Zero,BUF_FT(1,1,J),NP)
!end do
call DGEMM_('T','T',NP,NW*NJ,NB,One,CMOBLK,NB,BUF_HT,NW*NJ,Zero,BUF_FT,NP)

end subroutine FULLTRNSF
