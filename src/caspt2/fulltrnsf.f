************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE FULLTRNSF(NP,NW,NB,CMOBLK,NJ,BUF_HT,BUF_FT)
      IMPLICIT NONE
      INTEGER NP,NW,NB,NJ
      REAL*8 CMOBLK(NB,NP)
* In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)
      REAL*8 BUF_HT(NW*NJ,NB)
      REAL*8 BUF_FT(NP,NW*NJ)
*      DO J=1,NJ
** Compute fully transformed Cholesky vector buffer:
**  FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
*        CALL DGEMM_('T','T',NP,NW,NB,1.0D0,CMOBLK,NB,
*     &       BUF_HT(1,J,1),NW*NJ,0.0D0,BUF_FT(1,1,J),NP)
*      END DO
        CALL DGEMM_('T','T',NP,NW*NJ,NB,1.0D0,CMOBLK,NB,
     &       BUF_HT,NW*NJ,0.0D0,BUF_FT,NP)
      RETURN
      END

C =SVC= Special FULLTRNSF routine for boxed ordering with less
C efficiency than the original (except when we would take NWSZ=NW, which
C would make it possible to treat all J in one go)
      SUBROUTINE FULLTRNSF_BOXED(IPSTA,IWSTA,NPSZ,NWSZ,NP,NW,NB,
     &                           CMOBLK,NJ,BUF_HT,BUF_FT,BUFFY)
      IMPLICIT NONE
      INTEGER IPSTA,IWSTA,NPSZ,NWSZ,NP,NW,NB,NJ
      REAL*8 CMOBLK(NB,NPSZ)
C In: Vectors of type HALF(W,J,B) = Sum(CHO(AB,J)*CMO(A,W),A=1,NBAS)
      REAL*8 BUF_HT(NW,NJ,NB)
      REAL*8 BUF_FT(NP*NW,NJ)
      REAL*8 BUFFY(NPSZ,NWSZ)
      INTEGER IPW,J
C Compute fully transformed Cholesky vector buffer:
C FULL(P,W,J)=Sum(CMO(B,P)*HALF(W,J,B),B=1,NB)
      iPW=1+NW*(IPSTA-1)+NPSZ*(IWSTA-1)
      DO J=1,NJ
       CALL DGEMM_('T','T',NPSZ,NWSZ,NB,1.0D0,CMOBLK,NB,
     &      BUF_HT(IWSTA,J,1),NW*NJ,0.0D0,BUFFY,NPSZ)
       CALL DCOPY_(NPSZ*NWSZ,BUFFY,1,BUF_FT(iPW,J),1)
      ENDDO
      RETURN
      END
