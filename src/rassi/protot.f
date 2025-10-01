************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1999, Per Ake Malmqvist                                *
************************************************************************
      SUBROUTINE PROTOT(NPORB,NPSDSZ,IPSDMS,NPCSFSZ,IPCSFCP,PCSFTOSD)
      use rassi_aux, only: ipglob
      IMPLICIT REAL*8 (A-H,O-Z)
      Integer IPSDMS(NPORB,NPSDSZ),IPCSFCP(NPORB,NPCSFSZ)
      REal*8 PCSFTOSD(NPSDSZ,NPCSFSZ)
      INTEGER ASPIN, BSPIN, UPCPL, DWNCPL
      PARAMETER (UPCPL=1, DWNCPL=0)
      PARAMETER (ASPIN=1,BSPIN=0)
* Expand csf's in terms of determinants by the Grabenstetter method
*  (I.J.Q.C. 10, P142 (1976))
* Recoded by PAM 1999, after Jeppe Olsen.
*
* Input :
*         NPORB    : NUMBER OF OPEN ORBITALS
*         IPSDMS  : OCCUPATION OF PROTO-SDs
*         NPSDSZ   : NUMBER OF PROTOSD''s
*         NPCSFSZ  : NUMBER OF PROTOCSF''s
*         IPCSFCP  : SPIN COUPLINGS IN PROTO-CSFs
* Output :
*         PCSFTOSD :  NPSDSZ X NPCSFSZ MATRIX
*                GIVING EXPANSION FROM P-SD'S TO P-CSF'S


      DO JCSF = 1, NPCSFSZ
       IF( IPGLOB .GE. 5 ) WRITE(6,*) ' ....Output for P-CSF ',JCSF
       DO JDET = 1, NPSDSZ
C EXPANSION COEFFICIENT OF DETERMINANT JDET FOR P-CSF JCSF
        COEF1=1.0D0
        COEF2=1.0D0
        INDSMM=0
        INDSPM=0
        DO IOPEN = 1, NPORB
         IC=0
         IF(IPCSFCP(IOPEN,JCSF).EQ.UPCPL) IC=1
         IM=0
         IF(IPSDMS(IOPEN,JDET).EQ.ASPIN) IM=1
         IF(IC.EQ.0) THEN
          IF(IM.EQ.0) THEN
            INDSPM=INDSPM-1
            COEF1=COEF1*SQRT(DBLE(INDSPM+1))
* If COEF1 has gone down to 0 exactly.
            IF (INDSPM+1.eq.0) exit
          ELSE
            INDSMM=INDSMM-1
            COEF1=-COEF1*SQRT(DBLE(INDSMM+1))
* If COEF1 has gone down to 0 exactly.
            IF (INDSMM+1.eq.0) exit
          END IF
          COEF2=COEF2*SQRT(DBLE(INDSPM+INDSMM+2))
         ELSE
          IF(IM.EQ.0) THEN
            INDSMM=INDSMM+1
            COEF1=COEF1*SQRT(DBLE(INDSMM))
          ELSE
            INDSPM=INDSPM+1
            COEF1=COEF1*SQRT(DBLE(INDSPM))
          END IF
          COEF2=COEF2*SQRT(DBLE(INDSPM+INDSMM))
         END IF
        END DO

        PCSFTOSD(JDET,JCSF)=COEF1/COEF2

       END DO
      END DO

      END SUBROUTINE PROTOT
