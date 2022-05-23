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
#include "strbas.fh"
#include "csm.fh"
#include "stinf.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "WrkSpc.fh"
*
* ====================
*. Output common block : XISPSM is calculated
* ====================
*
#include "cicisp.fh"
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
      CALL GETMEM('KLBLTP','ALLO','INTE',KLBLTP,NSMST)
      KLCVST=1
c      IF(IDC.EQ.3 .OR. IDC .EQ. 4 )
c     &CALL MEMMAN(KLCVST,NSMST,'ADDL  ',2,'KLCVST')
      CALL GETMEM('KLIOIO','ALLO','INTE',KLIOIO,NOCTPA*NOCTPB)
*. Obtain array giving symmetry of sigma v reflection times string
*. symmetry.
c      IF(IDC.EQ.3.OR.IDC.EQ.4)
c     &CALL SIGVST(WORK(KLCVST),NSMST)

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
      CALL IAIBCM(ICI,iWORK(KLIOIO))

      DO  50 ISYM = 1, NSMCI
          CALL ZBLTP(ISMOST(1,ISYM),NSMST,IDC,
     &               iWORK(KLBLTP),iWORK(KLCVST))
          CALL NGASDT(IGSOCCX(1,1,ICI),IGSOCCX(1,2,ICI),
     &                NGAS,ISYM,NSMST,NOCTPA,NOCTPB,
     &                iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &                ISPGPFTP(1,IBSPGPFTP(IATP)),
     &                ISPGPFTP(1,IBSPGPFTP(IBTP)),
     &                MXPNGAS,NCOMB,XNCOMB,MXS,MXSOO,
     &                iWORK(KLBLTP),NTTSBL,LCOL,
     &                iWORK(KLIOIO),MXSOO_AS)
*

          XISPSM(ISYM,ICI) = XNCOMB
          MXSOOB = MAX(MXSOOB,MXSOO)
          MXSB = MAX(MXSB,MXS)
          MXSOOB_AS = MAX(MXSOO_AS,MXSOOB_AS)
          NBLKIC(ISYM,ICI) = NTTSBL
          LCOLIC(ISYM,ICI) = LCOL
   50 CONTINUE
      CALL GETMEM('KLBLTP','FREE','INTE',KLBLTP,NSMST)
      CALL GETMEM('KLIOIO','FREE','INTE',KLIOIO,NOCTPA*NOCTPB)
  100 CONTINUE
*
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
            WRITE(6,'(1H , 4E22.15)') (XISPSM(II,ICI),II=1,NSMCI)
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
*. Largest number of BLOCKS in a CI expansion
      MXNTTS = 0
      DO ICI = 1,NCMBSPC
       DO ISM =1, NSMCI
        MXNTTS = MAX(MXNTTS,NBLKIC(ISM,ICI))
       END DO
      END DO
*
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
*
*
      RETURN
      END
