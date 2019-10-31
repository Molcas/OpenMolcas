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
* Copyright (C) 1995,1999, Jeppe Olsen                                 *
************************************************************************
      SUBROUTINE PART_CIV2(     IDC,   IBLTP,   NSSOA,   NSSOB,  NOCTPA,
     &                       NOCTPB,   NSMST,   MXLNG,   IOCOC,  ISMOST,
     &                       NBATCH,  LBATCH, LEBATCH, I1BATCH,  IBATCH,
     &                        ICOMP, ISIMSYM)
*
* Jeppe Olsen
*
* Last update : May 1999 : ISIMSYM added
*
* Partition a CI vector into batches of blocks.
* The length of a batch must be at most MXLNG
* If ISIMSYM.EQ.1, TTS blocks that differs only in symmetry are not
*                  split.
*
* IF ICOMP. eq. 1 the complete civector is constructed
*
*
* Compared to PART_CIV, the NOOS arrays have been eliminated.
* They are becoming the size defining arrays - atleast at
* the laptop
*
*. Output
* NBATCH : Number of batches
* LBATCH : Number of blocks in a given batch
* LEBATCH : Number of elements in a given batch ( packed ) !
* I1BATCH : Number of first block in a given batch
* IBATCH : TTS blocks in Start of a given TTS block with respect to start
*          of batch
*   IBATCH(1,*) : Alpha type
*   IBATCH(2,*) : Beta sym
*   IBATCH(3,*) : Sym of alpha
*   IBATCH(4,*) : Sym of beta
*   IBATCH(5,*) : Offset of block with respect to start of block in
*                 expanded form
*   IBATCH(6,*) : Offset of block with respect to start of block in
*                 packed form
*   IBATCH(7,*) : Length of block, expandend form
*   IBATCH(8,*) : Length of block, packed form
*
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*.Input
C     INTEGER NOOS(NOCTPA,NOCTPB,NSMST)
C     INTEGER NOOSP(NOCTPA,NOCTPB,NSMST)
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
*.Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
C Dummy initialize
      INCLUDE = 0
      LBLOCKP = 0
*
      NTEST = 00
      IF(NTEST.GE.100) THEN
         WRITE(6,*)
         WRITE(6,*) ' =================='
         WRITE(6,*) '     PART_CIV2     '
         WRITE(6,*) ' =================='
         WRITE(6,*) ' IDC = ', IDC
         WRITE(6,*)
         WRITE(6,*) ' IOCOC Array '
         CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
         WRITE(6,*) ' ISMOST array '
         CALL IWRTMA(ISMOST,1,NSMST,1,NSMST)
         WRITE(6,*) ' IBLTP array '
         CALL IWRTMA(IBLTP,1,NSMST,1,NSMST)
         WRITE(6,*) ' NSSOA, NSSOB '
         CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
         CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
         WRITE(6,*) 'ISIMSYM, ICOMP = ', ISIMSYM,ICOMP
      END IF
*
*. block 1
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
*. New order : ISM,IB,IA (leftmost inner loop )
         IF(ISM.LT.NSMST) THEN
            ISM = ISM + 1
         ELSE
            ISM = 1
            IF(IB.LT.NOCTPB) THEN
               IB = IB + 1
            ELSE
               IB = 1
               IF(IA.LT.NOCTPA) THEN
                  IA = IA + 1
               ELSE
                  IFINI = 1
               END IF
            END IF
         END IF
      END IF
      IFRST = 0
      IF(IFINI.EQ.1) GOTO 2002
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. Size of TT block ( all symmetries)
      LBLOCK_AS = 0
c      IF(ISIMSYM.EQ.1 .AND. ISM. EQ. 1 ) THEN
c         DO IASM = 1, NSMST
c            IBSM = ISMOST(IASM)
c            NSTA = NSSOA(IASM,IA)
c            NSTB = NSSOB(IBSM,IB)
c            IF(IBLTP(IASM).EQ.0) GOTO 99
c            IF(IBLTP(IASM).EQ.2.AND.IA.LT.IB) GOTO 99
c            LBLOCK_AS = LBLOCK_AS + NSTA*NSTB
c 99         CONTINUE
c         END DO
c         INCLUDE = 0
C?      WRITE(6,*) ' IA IB LBLOCK_AS', IA,IB, LBLOCK_AS
c         IF(LENGTH+LBLOCK_AS.LE.MXLNG.OR.ICOMP.EQ.1) INCLUDE = 1
c      END IF
*. Should this block be included
      IASM = ISM
      IBSM = ISMOST(IASM)
      IF(IDC.EQ.2) THEN
         IF(IA.LT.IB) GOTO 1000
         IF(IA.EQ.IB.AND.IASM.LT.IBSM) GOTO 1000
      END IF
*. can this block be included
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
      IF(IDC.EQ.1.OR.IA.GT.IB.OR.(IA.EQ.IB.AND.IASM.GT.IBSM)) THEN
         LBLOCKP = NSTA*NSTB
      ELSE IF(IDC.EQ.2.AND.IA.EQ.IB.AND.IASM.EQ.IBSM) THEN
         LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C?      WRITE(6,*) ' IASM IBSM IA IB LBLOCKP,LBLOCK' ,
C?   &               IASM,IBSM,IA,IB,LBLOCKP,LBLOCK
*
c      IF(ISIMSYM.EQ.0) THEN
         INCLUDE = 0
         IF(LENGTH+LBLOCK.LE.LBLOCK.OR.ICOMP.EQ.1) INCLUDE = 1
c      END IF
*
      IF(INCLUDE.EQ.1) THEN
         NBLOCK = NBLOCK + 1
         IBLOCK = IBLOCK + 1
         LBATCH(NBATCH) = LBATCH(NBATCH)+1
         IBATCH(1,IBLOCK) = IA
         IBATCH(2,IBLOCK) = IB
         IBATCH(3,IBLOCK) = ISM
         IBATCH(4,IBLOCK) = IBSM
         IBATCH(5,IBLOCK) = LENGTH+1
         IBATCH(6,IBLOCK) = LENGTHP+1
         IBATCH(7,IBLOCK) = LBLOCK
         IBATCH(8,IBLOCK) = LBLOCKP
         LENGTH = LENGTH + LBLOCK
         LENGTHP= LENGTHP+ LBLOCKP
         LEBATCH(NBATCH) = LENGTHP
         GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.INCLUDE.EQ.0.AND.NBLOCK.EQ.0) THEN
         WRITE(6,*)
     &        ' Not enough space to include a single Block'
         WRITE(6,*) ' Since I cannot procede I will stop '
         WRITE(6,*) ' Insufficient space detected in PART_CIV'
         WRITE(6,*) ' Alter GAS space or raise Buffer from ', MXLNG
*         STOP 'Error in PART_CIV2'
         CALL SYSABENDMSG('lucia_util/part_civ2','Internal error',' ')
      ELSE
*. This batch is finished, goto next batch
         GOTO 2000
      END IF
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
C?      WRITE(6,*) 'Output from PART_CIV'
C?      WRITE(6,*) '====================='
         WRITE(6,*)
         WRITE(6,*) ' Number of batches ', NBATCH
         DO JBATCH = 1, NBATCH
            WRITE(6,*)
            WRITE(6,*) ' Info on batch ', JBATCH
            WRITE(6,*) ' *********************** '
            WRITE(6,*)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lengths of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
               WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
            END DO
         END DO
      END IF
*
      RETURN
      END
