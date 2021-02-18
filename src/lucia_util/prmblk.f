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
* Copyright (C) 1993, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE PRMBLK(     IDC,    ISGV,    IASM,    IBSM,    IATP,
     &                      IBTP,      PS,      PL,    JATP,    JBTP,
     &                      JASM,    JBSM,    ISGN,    ITRP,   NPERM)
*
* A block of CI coefficients defined by by IATP,IASM,IBTP,IBSM is given
*
* Obtain the number of other blocks that can be obtained by spin
* and relection symmetry.
*
* Jeppe Olsen, July 1993
*
* =====
* Output
* =====
* JATP(I),JASM(I),JBTP(I),JBSM(I) Indeces for Block I
* NPERM : Number of blocks  that can be obtained
* ITRP(I) = 1 => block should     be transposed
*         = 0 => block should not be transposed
* ISGN   : Sign to multiply previous block with to getnew sign
*
*
* There are four types of permutations
*
*    operation   *      JASM  *      JBSM  * JATP * JBTP * Iperm * Sign *
*   *********************************************************************
*   * Identity   *      IASM  *      IBSM  * IATP * IBTP *   0   * 1    *
*   * Ml         * ISGV(IASM) * ISGV(IBSM) * IATP * IBTP *   0   * PL   *
*   * Ms         *      IBSM  *      IASM  * IBTP * IATP *   1   * PS   *
*   * Ms+Ml      * ISGV(IBSM) * ISGV(IASM) * IBTP * IATP *   1   * PS PL*
*   *********************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
*.Input
      DIMENSION ISGV(*)
*.Output
      DIMENSION JATP(4),JBTP(4),JASM(4),JBSM(4),ISGN(4),ITRP(4)
*
*. To eliminate some compiler warnings
      KASM = 0
      KBSM = 0
      KATP = 0
      KBTP = 0
      KSIGN = 0
      KTRP = 0
      LSIGN = 0
      LTRP = 0
*
      NPERM = 0
      DO 100 IPERM = 1, 4
        ISET = 0
        IF(IPERM.EQ.1) THEN
*
* Identity operation
*
          KASM = IASM
          KBSM = IBSM
          KATP = IATP
          KBTP = IBTP
          KSIGN = 1
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.2.AND.(IDC.EQ.3.OR.IDC.EQ.4)) THEN
*
* Ml reflection
*
          KASM = ISGV(IASM)
          KBSM = ISGV(IBSM)
          KATP = IATP
          KBTP = IBTP
          IF(PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PL .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 0
          ISET = 1
        ELSE IF(IPERM.EQ.3.AND.(IDC.EQ.2.OR.IDC.EQ.4)) THEN
*
* Ms reflection
*
          KASM = IBSM
          KBSM = IASM
          KATP = IBTP
          KBTP = IATP
          IF(PS.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        ELSE IF(IPERM.EQ.4 .AND. IDC.EQ.4) THEN
*
* Ms Ml  reflection
*
          KASM = ISGV(IBSM)
          KBSM = ISGV(IASM)
          KATP = IBTP
          KBTP = IATP
          IF(PS*PL.EQ.1.0D0) THEN
            KSIGN = 1
          ELSE IF (PS .EQ. -1.0D0) THEN
            KSIGN = -1
          END IF
          KTRP = 1
          ISET = 1
        END IF
*
        IF(ISET.EQ.1) THEN
*. A new permutation was found, check and see if it was obtained previously
          INEW = 1
          DO 50 LPERM = 1, NPERM
            IF(JATP(LPERM).EQ.KATP  .AND. JASM(LPERM).EQ.KASM .AND.
     &         JBTP(LPERM).EQ.KBTP  .AND. JBSM(LPERM).EQ.KBSM) INEW = 0
   50     CONTINUE
          IF(INEW.EQ.1) THEN
*. The permutation was new, add it to the list
            NPERM = NPERM + 1
            JASM(NPERM) = KASM
            JBSM(NPERM) = KBSM
            JATP(NPERM) = KATP
            JBTP(NPERM) = KBTP
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KSIGN.EQ.LSIGN))THEN
              ISGN(NPERM) = 1
            ELSE
              ISGN(NPERM) = -1
            END IF
            LSIGN = KSIGN
            IF(NPERM.EQ.1. OR. (NPERM.GE.1.AND.KTRP.EQ.LTRP))THEN
              ITRP(NPERM) = 0
            ELSE
              ITRP(NPERM) = 1
            END IF
            LTRP = KTRP
          END IF
        END IF
  100 CONTINUE
*
*. Should the block be trnasposed or scaled to return to initial form
      ITRP(NPERM+1) = LTRP
      ISGN(NPERM+1) = LSIGN
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,'(A,4I4)') ' Blocks obtained from IASM IBSM IATP IBTP ',
     &  IASM,IBSM,IATP,IBTP
        WRITE(6,*)
        WRITE(6,'(A)') ' JASM JBSM JATP JBTP Isgn Itrp  '
        WRITE(6,*)
        DO 10 IPERM = 1, NPERM
          WRITE(6,'(2x,6I4)') JASM(IPERM),JBSM(IPERM),JATP(IPERM),
     &                        JBTP(IPERM),ISGN(IPERM),ITRP(IPERM)
   10   CONTINUE
      END IF
*
      RETURN
      END
