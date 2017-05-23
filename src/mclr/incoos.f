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
* Copyright (C) 1991, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE INCOOS(IDC,IBLTP,NOOS,NOCTPA,NOCTPB,ISTSM,ISTTA,ISTTB,
     &                  NSMST,IENSM,IENTA,IENTB,IACOOS,MXLNG,IFINI,
     &                  NBLOCK,INCFST,IOCOC)
*
* Obtain Number of OOS blocks that can be included
* IN MXLNG word starting from block after ISTSM,ISTTA,ISTTB
* Activated blocks are given in IACOOS
* Last activated block is (IENSM,IENTA,IENTB)
* If all blocks have been accessed IFINI is returned as 1
* Diagonal blocks are expanded
*
* Jeppe Olsen, Winter of 1991
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
      INTEGER NOOS(NOCTPA,NOCTPB,NSMST)
      INTEGER IOCOC(NOCTPA,NOCTPB)
C-May 7
      INTEGER IBLTP(*)
C-May 7
*.Output
      INTEGER IACOOS(NOCTPA,NOCTPB,NSMST)
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' =================='
        WRITE(6,*) ' INCOOS in action  '
        WRITE(6,*) ' =================='
        WRITE(6,*)
        WRITE(6,*) ' NOOS(NOCTPA,NOCTPB,NSMST) array (input) '
        WRITE(6,*)
        DO ISMST = 1, NSMST
         WRITE(6,*) ' ISMST = ', ISMST
         CALL IWRTMA(NOOS(1,1,ISMST),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        END DO
      END IF
*
      IPA = 0
      IPB = 0
      IPSM = 0
*
*.Initialize
      CALL ISETVC(IACOOS,0,NOCTPA*NOCTPB*NSMST)
      IFRST = 1
      ISM = ISTSM
      IA = ISTTA
      IB = ISTTB
      LENGTH = 0
      NBLOCK = 0
      IENSM = ISTSM
      IENTA = ISTTA
      IENTB = ISTTB
      IFINI = 0
      IF(INCFST.EQ.1) GOTO 999
 1000 CONTINUE
*.Next block
      IPA = IA
      IPB = IB
      IPSM = ISM
*
      IF(IB.LT.NOCTPB) THEN
        IB = IB + 1
      ELSE
        IB = 1
        IF(IA.LT.NOCTPA) THEN
          IA = IA+ 1
        ELSE
          IA = 1
          IF(ISM.LT.NSMST) THEN
            ISM = ISM + 1
          ELSE
            IFINI = 1
          END IF
        END IF
      END IF
      IF(IFINI.EQ.1) GOTO 1001
*. Should this block be included
  999 CONTINUE
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.0) GOTO 1000
      IF(IDC.NE.1.AND.IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
C?    write(6,*) ' INCOOS IDC IBLTP ', IDC,IBLTP(ISM)
*. can this block be included
      LBLOCK = NOOS(IA,IB,ISM)
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG) THEN
        NBLOCK = NBLOCK + 1
        LENGTH = LENGTH + LBLOCK
        IACOOS(IA,IB,ISM) = 1
        IF(NBLOCK.EQ.1) THEN
          ISTTA = IA
          ISTTB = IB
          ISTSM = ISM
         END IF
        GOTO 1000
      ELSE
        IA = IPA
        IB = IPB
        ISM = IPSM
      END IF
 1001 CONTINUE
*
      IENSM = ISM
      IENTA = IA
      IENTB = IB
      IF(IFINI.EQ.0.AND.NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in INCOOS '
        WRITE(6,*) ' Alter RAS space of raise Buffer from ', MXLNG
        CALL MEMCHK
*        STOP 11
        CALL SYSABENDMSG('lucia_util/incoos','Internal error',' ')
      END IF
*
      IF(NTEST.NE.0) THEN
        WRITE(6,*) 'Output from INCOOS '
        WRITE(6,*) '==================='
        WRITE(6,*)
     &  ' Length and number of included blocks ',LENGTH,NBLOCK
      END IF
      IF(NTEST.GE.2) THEN
        DO 100 ISM = ISTSM,IENSM
          WRITE(6,*) ' Active blocks of symmetry ',ISM
          CALL IWRTMA(IACOOS(1,1,ISM),NOCTPA,NOCTPB,NOCTPA,NOCTPB)
  100   CONTINUE
        IF(IFINI.EQ.1) WRITE(6,*) ' No new blocks '
      END IF
*
      RETURN
      END
