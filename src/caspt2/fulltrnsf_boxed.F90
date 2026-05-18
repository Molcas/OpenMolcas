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
! =SVC= Special FULLTRNSF routine for boxed ordering with less
! efficiency than the original (except when we would take NWSZ=NW, which
! would make it possible to treat all J in one go)

subroutine FULLTRNSF_BOXED(IPSTA,IWSTA,NPSZ,NWSZ,NP,NW,NB,CMOBLK,NJ,BUF_HT,BUF_FT,BUFFY)

use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IPSTA, IWSTA, NPSZ, NWSZ, NP, NW, NB, NJ
real(kind=wp), intent(in) :: CMOBLK(NB,NPSZ), BUF_HT(NW,NJ,NB)
real(kind=wp), intent(out) :: BUF_FT(NP*NW,NJ), BUFFY(NPSZ,NWSZ)
integer(kind=iwp) :: IPW, J

! In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)

! Compute fully transformed Cholesky vector buffer:
! FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
iPW = 1+NW*(IPSTA-1)+NPSZ*(IWSTA-1)
do J=1,NJ
  call DGEMM_('T','T',NPSZ,NWSZ,NB,One,CMOBLK,NB,BUF_HT(IWSTA,J,1),NW*NJ,Zero,BUFFY,NPSZ)
  call DCOPY_(NPSZ*NWSZ,BUFFY,1,BUF_FT(iPW,J),1)
end do

end subroutine FULLTRNSF_BOXED
