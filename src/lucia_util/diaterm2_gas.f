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
* Copyright (C) 1995, Jeppe Olsen                                      *
************************************************************************
      SUBROUTINE DIATERM2_GAS( FACTOR,  ITASK,    VEC, NBLOCK, IBLOCK,
     &                           IOFF,  JPERT,    J12,    JDC)
* = DIATERM_GAS, just J12 added !
*
* Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
* Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
*
* For the NBLOCKS givem in IBLOCK starting from BLOCK IOFF
*
* If JPERT.NE.0, the perturbation operator as defined by IPART is used.
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.fh"
#include "orbinp.fh"
#include "cicisp.fh"
#include "strbas.fh"
#include "cstate.fh"
#include "strinp.fh"
#include "stinf.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "cprnt.fh"
#include "cgas.fh"
#include "gasstr.fh"
#include "glbbas.fh"
#include "oper.fh"
#include "cecore.fh"
*
#include "cintfo.fh"
*
      INTEGER IBLOCK(8,*)
*
      DIMENSION VEC(*)
*
      CALL QENTER('DIATR')
      IDUM=0
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP = 1
      IBTP = 2
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
C     IF(JPERT.EQ.0) THEN
*. Use full Hamiltonian
C       I12 = 2
C       IPERTOP = 0
C     ELSE
*. Use perturbation operator
C       IF(IPART.EQ.1) THEN
*. Moller-Plesset partitioning
C         I12 = 1
C         IPERTOP = 1
C       ELSE IF(IPART.EQ.2) THEN
*. Epstein-Nesbet Partitioning
C         I12 = 2
C         IPERTOP = 0
C       END IF
C     END IF

      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   DIATERM2_GAS speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        WRITE(6,*) ' JPERT,IPART,J12,IPERTOP',JPERT,J12,IPERTOP
      END IF
*. A bit of scracth
      CALL GETMEM('KLJ   ','ALLO','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','ALLO','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','ALLO','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXB  ','ALLO','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLH1D ','ALLO','REAL',KLH1D ,NACOB)
*. Space for blocks of strings
      CALL GETMEM('KLASTR','ALLO','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','ALLO','INTE',KLBSTR,MXNSTR*NAEL)
      MAXA = IMNMX(IWORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL GETMEM('KLRJKA','ALLO','REAL',KLRJKA,MAXA)
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      CALL GT1DIA(WORK(KLH1D))
C!    IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KFI),WORK(KINT1),NINT1)
      IF(J12.EQ.2)
     &CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS,IREOST)
*. Core energy not included
      ECOREP = 0.0D0
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
      CALL DIATERMS_GAS(NAEL,IWORK(KLASTR),NBEL,IWORK(KLBSTR),
     &                  NACOB,VEC,NSMST,
     &                  WORK(KLH1D),JDC,WORK(KLXB),WORK(KLJ),WORK(KLK),
     &                  iWORK(KNSTSO(IATP)),iWORK(KNSTSO(IBTP)),
     &                  ECOREP,0,0,IPRDIA,NTOOB,WORK(KLRJKA),J12,
     &                  IBLOCK(1,IOFF),NBLOCK,ITASK, FACTORX,0,[0])
*
C    &                  IBLOCK,NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
*.Flush local memory
      CALL GETMEM('KLJ   ','FREE','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','FREE','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','FREE','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXB  ','FREE','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLH1D ','FREE','REAL',KLH1D ,NACOB)
      CALL GETMEM('KLASTR','FREE','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','FREE','INTE',KLBSTR,MXNSTR*NAEL)
      CALL GETMEM('KLRJKA','FREE','REAL',KLRJKA,MAXA)
*
      CALL QEXIT('DIATR')
*
      IF(NTEST.GE.100) THEN
        WRITE(6,*)  ' output vector from DIATRM '
        CALL WRTTTS(      VEC,IBLOCK(1,IOFF),NBLOCK,  NSMST,
     &                 IWORK(KNSTSO(IATP)),IWORK(KNSTSO(IBTP)),IDC)
      END IF
*
      RETURN
      END
