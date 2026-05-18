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
      SUBROUTINE FULLTRNSF(NP,NW,NB,CMOBLK,NJ,BUF_HT,BUF_FT)
      use constants, only: Zero, One
      use definitions, only: iwp, wp
      IMPLICIT NONE
      INTEGER(kind=iwp), intent(in):: NP,NW,NB,NJ
      REAL(kind=wp), intent(in):: CMOBLK(NB,NP)
! In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)
      REAL(kind=wp), intent(in):: BUF_HT(NW*NJ,NB)
      REAL(kind=wp), intent(out):: BUF_FT(NP,NW*NJ)
!      DO J=1,NJ
!* Compute fully transformed Cholesky vector buffer:
!*  FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
!        CALL DGEMM_('T','T',NP,NW,NB,One,CMOBLK,NB,
!     &       BUF_HT(1,J,1),NW*NJ,Zero,BUF_FT(1,1,J),NP)
!      END DO
        CALL DGEMM_('T','T',NP,NW*NJ,NB,One,CMOBLK,NB,                  &
     &       BUF_HT,NW*NJ,Zero,BUF_FT,NP)
      END SUBROUTINE FULLTRNSF
