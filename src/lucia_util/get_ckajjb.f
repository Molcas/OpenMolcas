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
      SUBROUTINE GET_CKAJJB(     CB,     NJ,    NJA, CKAJJB,    NKA,
     &                          NJB,      J,   ISCA,   SSCA)
*
* Obtain for given orbital index j the gathered matrix
*
* C(Ka,j,Jb) = SSCA(Ka)C(Jb,ISCA(Ka))
*
* For efficient processing of alpha-beta loop
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
       DIMENSION CB(NJB,NJA), SSCA(*),ISCA(*)
*. Output
       DIMENSION CKAJJB(*)
*
*. To get rid of annoying and incorrect compiler warnings
      ICOFF = 0
*
C?    WRITE(6,*) ' From GET_CKAJJB'
C     LBLK = 100
      LBLK = 40
      NBLK = NJB/LBLK
      IF(LBLK*NBLK.LT.NJB) NBLK = NBLK + 1
      DO ICBL = 1, NBLK
        IF(ICBL.EQ.1) THEN
          ICOFF = 1
        ELSE
          ICOFF = ICOFF + LBLK
        END IF
        ICEND = MIN(ICOFF+LBLK-1,NJB)
        ICONST = NKA*NJ
        IADR0 =  (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
        IF(ICEND.GT.ICOFF) THEN
*. Inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
                CKAJJB(IADR) = S*CB(JB,IROW)
              END DO
            ELSE
              IADR = IADR0 + KA
              DO JB = ICOFF,ICEND
C               IADR = KA + (J-1)*NKA+(JB-1)*NKA*NJ
                IADR = IADR + ICONST
                CKAJJB(IADR) = 0.0D0
              END DO
            END IF
          END DO
        ELSE
*. No inner loop over JB
          DO KA  = 1, NKA
            IF(ISCA(KA).NE.0) THEN
              S = SSCA(KA)
              IROW = ISCA(KA)
C             IADR = KA + (J-1)*NKA+(ICOFF-1-1)*NKA*NJ
              IADR = IADR0 + KA
C             DO JB = ICOFF,ICEND
*. Adress of C(Ka,j,Jb)
                IADR = IADR + ICONST
                CKAJJB(IADR) = S*CB(ICOFF,IROW)
C             END DO
            ELSE
              IADR = IADR0 + KA
C             DO JB = ICOFF,ICEND
C               IADR = KA + (J-1)*NKA+(JB-1)*NKA*NJ
                IADR = IADR + ICONST
                CKAJJB(IADR) = 0.0D0
C             END DO
            END IF
          END DO
        END IF
*       ^ End of test ICEND,ICOFF
      END DO
*
      RETURN
      END
*
