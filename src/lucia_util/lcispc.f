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
* Copyright (C) 1994,1995,1999, Jeppe Olsen                            *
************************************************************************
      SUBROUTINE LCISPC(IPRNT)
      use stdalloc, only: mma_allocate, mma_deallocate
      use strbas
*
* Number of dets and combinations
* per symmetry for each type of internal space
*
* Jeppe Olsen , Winter 1994/1995 ( woops !)
*               MXSOOB_AS added,MXSB removed May 1999
*
* GAS VERSION
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* ===================
*.Input common blocks
* ===================
*
#include "mxpdim.fh"
#include "lucinp.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "csm.fh"
#include "stinf.fh"
#include "cgas.fh"
#include "gasstr.fh"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.fh"

      Integer, Allocatable:: LBLTP(:), LIOIO(:), CVST(:)
*
*
*. Number of spaces
      NICISP = NCMBSPC
C?    write(6,*) ' LCISPC : NICISP ', NICISP
*. Type of alpha- and beta strings
      IATP = 1
      IBTP = 2
*
      NOCTPA =  NOCTYP(IATP)
      NOCTPB =  NOCTYP(IBTP)
*.Local memory
      CALL mma_allocate(LBLTP,NSMST,Label='LBLTP')
      Call mma_allocate(CVST,NSMST,Label='CVST')
      CALL mma_allocate(LIOIO,NOCTPA*NOCTPB,Label='LIOIO')
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
c      IF(IDC.EQ.3.OR.IDC.EQ.4) CALL SIGVST(CVST,NSMST)

*. Array defining symmetry combinations of internal strings
*. Number of internal dets for each symmetry
        CALL SMOST(NSMST,NSMCI,MXPCSM,ISMOST)
*. MXSB is not calculated anymore, set to 0
      MXSB = 0
*
      MXSOOB = 0
      MXSOOB_AS = 0
      DO 100 ICI = 1, NICISP
*. allowed combination of types
      CALL IAIBCM(ICI,LIOIO)

      DO  50 ISYM = 1, NSMCI
          CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,LBLTP,CVST)
          CALL NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),
     &                NGAS,ISYM,NSMST,NOCTPA,NOCTPB,
     &                NSTSO(IATP)%I,NSTSO(IBTP)%I,
     &                ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                ISPGPFTP(1,IBSPGPFTP(IBTP)),
     &                MXPNGAS,NCOMB,XNCOMB,MXS,MXSOO,
     &                LBLTP,NTTSBL,LCOL,
     &                LIOIO,MXSOO_AS)
*

          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB = MAX(MXSOOB,MXSOO)
          MXSB = MAX(MXSB,MXS)
          MXSOOB_AS = MAX(MXSOO_AS,MXSOOB_AS)
          NBLKIC(ISYM,ICI) = NTTSBL
          LCOLIC(ISYM,ICI) = LCOL
   50 CONTINUE
      Call mma_deallocate(LBLTP)
      Call mma_deallocate(CVST)
      Call mma_deallocate(LIOIO)
  100 CONTINUE
*
#ifdef _DEBUGPRINT_
      NTEST = 0
      NTEST = MAX(NTEST,IPRNT)
      IF (NTEST .GE. 5) THEN
         WRITE(6,*)
         WRITE(6,*)
         WRITE(6,*)
     &      ' Number of internal combinations per symmetry '
         WRITE(6,*)
     &      ' =========================================== '
*
         DO 200 ICI = 1, NCMBSPC
            WRITE(6,*) ' CI space ', ICI
            WRITE(6,'(1X, 4ES22.15)') (XISPSM(II,ICI),II=1,NSMCI)
C         CALL WRTMAT(XISPSM(1,ICI),1,NSMCI,1,NSMCI)
  200    CONTINUE
         WRITE(6,*)
         WRITE(6,*) ' Largest Symmetry-type-type block ',MXSOOB
         WRITE(6,*) ' Largest type-type block (all symmetries) ',
     &      MXSOOB_AS
         WRITE(6,*)
*
         WRITE(6,*)
     &      ' Number of TTS subblocks per CI expansion '
         WRITE(6,*)
     &      ' ======================================== '
*
        DO  ICI = 1,  NCMBSPC
            WRITE(6,*) ' Internal CI space ', ICI
            CALL IWRTMA(NBLKIC(1,ICI),1,NSMCI,1,NSMCI)
        END DO
      END IF
#else
      Call Unused_Integer(IPRNT)
#endif
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
#ifdef _DEBUGPRINT_
      IF(NTEST.GE.5) THEN
      WRITE(6,*) ' Largest number of blocks in CI expansion',
     &   MXNTTS
*
      WRITE(6,*)
     &' Number of columns per CI expansion '
      WRITE(6,*)
     & ' =================================== '
*
      DO  ICI = 1,  NCMBSPC
          WRITE(6,*) ' Internal CI space ', ICI
          CALL IWRTMA(LCOLIC(1,ICI),1,NSMCI,1,NSMCI)
      END DO
      END IF
#endif
*
      END
