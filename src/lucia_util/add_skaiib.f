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
      SUBROUTINE ADD_SKAIIB(     SB,     NI,    NIA, SKAIIB,    NKA,
     &                          NIB,      I,   ISCA,   SSCA)
*
* Update Transposed sigma block with contributions for given orbital index j
* from the matrix S(Ka,i,Ib)
*
* S(Ib,Isca(Ka)) =  S(Ib,Isca(Ka)) + Ssca(Ka)*S(Ka,I,Ib)
*
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION SKAIIB(*),SSCA(*),ISCA(*)
*. Input and Output
       DIMENSION SB(NIB,NIA)
*
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
*
C     LBLK = 100
      LBLK = 40
      NBLK = NIB/LBLK
      IF(LBLK*NBLK.LT.NIB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NIB)
        ICONST = NKA*NI
        IADR0 =  (I-1)*NKA+(ICOFF-1-1)*NKA*NI
        IF(ICEND.GT.ICOFF) THEN
*. Use form with Inner loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (I-1)*NKA+(ICOFF-1-1)*NKA*NI
              IADR = IADR0 + KA
              DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
                IADR = IADR + ICONST
                SB(Ib,IROW) = SB(Ib,IROW)+S*SKAIIB(IADR)
              END DO
            END IF
          END DO
        ELSE
*. Form with no loop over IB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
              IADR = IADR0 + KA + ICONST
C             DO IB = ICOFF,ICEND
*. Adress of S(Ka,i,Ib)
C               IADR = IADR + ICONST
                SB(ICOFF,IROW) = SB(ICOFF,IROW)+S*SKAIIB(IADR)
C             END DO
            END IF
          END DO
        END IF
*       ^ End of test of ICOFF=ICEND
      END DO
*
      RETURN
      END
C               GET_CKAJJB(CB,NJ,NJA,CJRES,NKABTC,NJB,
C    &                          JJ,I1(1,JJ),XI1S(1,JJ)
