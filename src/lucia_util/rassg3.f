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
* Copyright (C) 1991,1997, Jeppe Olsen                                 *
*               2015, Lasse Kragh Soerensen                            *
************************************************************************
      SUBROUTINE RASSG3(      CB,      SB,   NBATS,   LBATS,  LEBATS,
     &                    I1BATS,   IBATS,     LUC,    LUHC,
     &                    I_AM_OUT,N_ELIMINATED_BATCHES)
*
* Direct RAS routine employing combined MOC/n-1 resolution method
*
* Jeppe Olsen   Winter of 1991
*               May 1997 : Connected to SBLOCK
*
* Lasse Soerensen October 2015
*                 Do not calculate unwanted batches for highly
*                 excited states.
*
* =====
* Input
* =====
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "para_info.fh"
#include "WrkSpc.fh"
#include "io_util.fh"
*. Batches of sigma
      INTEGER LBATS(*),LEBATS(*),I1BATS(*),IBATS(8,*)
      DIMENSION I_AM_OUT(*)
*.Scratch
      DIMENSION SB(*),CB(*)
*
      CALL QENTER('RASSG')
      NTEST = 00
C     NTEST = MAX(NTEST,IPRNT)
      IF(NTEST.GE.20) THEN
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 speaking :'
        WRITE(6,*) ' ================='
        WRITE(6,*) ' RASSG3 : NBATS = ',NBATS
      END IF
*
CSVC: Compute offsets of a sigma batch in the sigma array.
C     The batches used inside sblock(s) use a batch size corresponding
C     to the 'expanded form' as computed inside part_civ2. This is
C     stored inside 7th element of IBATS. Later, the size that needs to
C     be actually written to disc uses the 'packed form', stored inside
C     the 8th element of IBATS. This also computes the total size NSB.
      CALL GETMEM('SBSIZ','ALLO','INTE',LSBSIZ,NBATS)
      CALL GETMEM('SBOFF','ALLO','INTE',LSBOFF,NBATS)
      NSB=0
      DO JBATS = 1, NBATS
        ISTA=I1BATS(JBATS)
        IEND=I1BATS(JBATS)+LBATS(JBATS)-1
        IWORK(LSBSIZ+JBATS-1)=SUM(IBATS(7,ISTA:IEND))
        NSB=NSB+IWORK(LSBSIZ+JBATS-1)
      END DO
      IWORK(LSBOFF) = 1
      DO JBATS = 2, NBATS
        IWORK(LSBOFF+JBATS-1) = IWORK(LSBOFF+JBATS-2) +
     &                           IWORK(LSBSIZ+JBATS-2)
      END DO

CSVC: the entire sigma array is zeroed here, because each process will
C     zero only its own sigma blocks, and we need to do a global sum
C     operations later to combine blocks before writing.
      CALL DCOPY_(NSB,0.0D0,0,SB,1)

      DO JBATS=1,NBATS
*
* Lasse addition start
*
        I_AM_NOT_WANTED = 0
        DO I = 1, N_ELIMINATED_BATCHES
          IF(I_AM_OUT(I).EQ.JBATS) THEN
            I_AM_NOT_WANTED = 1
            EXIT
          END IF
        END DO
        IF(I_AM_NOT_WANTED.EQ.1) CYCLE
*
* Lasse addition end
*

      ISBOFF=IWORK(LSBOFF+JBATS-1)
*. Obtain sigma for batch of blocks
      CALL SBLOCK(LBATS(JBATS),IBATS(1,I1BATS(JBATS)),1,
     &            CB,SB(ISBOFF),LUC,0,0,0,0,0)

      END DO

      CALL GADSUM(SB,NSB)
CSVC: Write sigma array to disk here, after sum reduction.
C     The writing is done in consecutive blocks, but since I don't know
C     if this block structure is used internally, I didn't optimize this.
      IF(LUHC.GT.0) IDISK(LUHC)=0
      DO JBATS = 1, NBATS
        ISBOFF=IWORK(LSBOFF+JBATS-1)
        DO ISBLK = I1BATS(JBATS),I1BATS(JBATS)+ LBATS(JBATS)-1
          IOFF = IBATS(6,ISBLK)
          ILEN = IBATS(8,ISBLK)
          CALL ITODS(ILEN,1,-1,LUHC)
          CALL TODSC(SB(ISBOFF-1+IOFF),ILEN,-1,LUHC)
        END DO
      END DO

      CALL GETMEM('SBOFF','FREE','INTE',LSBOFF,NBATS)
      CALL GETMEM('SBSIZ','FREE','INTE',LSBSIZ,NBATS)

      CALL ITODS(-1,1,-1,LUHC)
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Final S-vector on disc'
        CALL WRTVCD(SB,LUHC,1,-1)
      END IF
*
      CALL QEXIT('RASSG')
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) Call Unused_integer_array(LEBATS)
      END
