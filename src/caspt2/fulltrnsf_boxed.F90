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
      SUBROUTINE FULLTRNSF_BOXED(IPSTA,IWSTA,NPSZ,NWSZ,NP,NW,NB,        &
     &                           CMOBLK,NJ,BUF_HT,BUF_FT,BUFFY)
      use constants, only: Zero, One
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp) IPSTA,IWSTA,NPSZ,NWSZ,NP,NW,NB,NJ
      REAL(kind=wp), intent(in):: CMOBLK(NB,NPSZ)
! In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)
      REAL(kind=wp), intent(in):: BUF_HT(NW,NJ,NB)
      REAL(kind=wp), intent(out):: BUF_FT(NP*NW,NJ)
      REAL(kind=wp), intent(out):: BUFFY(NPSZ,NWSZ)
      INTEGER(kind=iwp) IPW,J
! Compute fully transformed Cholesky vector buffer:
! FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
      iPW=1+NW*(IPSTA-1)+NPSZ*(IWSTA-1)
      DO J=1,NJ
       CALL DGEMM_('T','T',NPSZ,NWSZ,NB,One,CMOBLK,NB,                  &
     &      BUF_HT(IWSTA,J,1),NW*NJ,Zero,BUFFY,NPSZ)
       CALL DCOPY_(NPSZ*NWSZ,BUFFY,1,BUF_FT(iPW,J),1)
      ENDDO
      END SUBROUTINE FULLTRNSF_BOXED
