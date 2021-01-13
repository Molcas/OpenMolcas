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
      SUBROUTINE GASDIAT(    DIAG,   LUDIA,   ECORE,  ICISTR,     I12,
     &                      IBLTP,  NBLOCK,  IBLKFO)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
      IMPLICIT REAL*8(A-H,O-Z)
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
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
#include "io_util.fh"
*
      DIMENSION IBLTP(*)
      DIMENSION IBLKFO(8,NBLOCK)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
*
      NTEST = 0
      NTEST = MAX(NTEST,IPRDIA)
*
** Specifications of internal space
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
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ================'
        WRITE(6,*) ' GASDIA speaking '
        WRITE(6,*) ' ================'
        WRITE(6,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
        write(6,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
        write(6,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
      END IF
*
**. Local memory
*
      CALL GETMEM('KLJ   ','ALLO','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','ALLO','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','ALLO','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXB  ','ALLO','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLH1D ','ALLO','REAL',KLH1D ,NACOB)
*. Space for blocks of strings
      CALL GETMEM('KLASTR','ALLO','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','ALLO','INTE',KLBSTR,MXNSTR*NBEL)
      MAXA = IMNMX(IWORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL GETMEM('KLRJKA','ALLO','REAL',KLRJKA,MAXA)
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA(WORK(KLH1D))
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS,IREOST)
      IF( LUDIA .GT. 0 ) IDISK(LUDIA)=0
      CALL GASDIAS(NAEL,IWORK(KLASTR),NBEL,IWORK(KLBSTR),
     &             NACOB,DIAG,NSMST,
     &             WORK(KLH1D),WORK(KLXB),WORK(KLJ),WORK(KLK),
     &             IWORK(KNSTSO(IATP)),IWORK(KNSTSO(IBTP)),
     &             LUDIA,ECORE,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &             WORK(KLRJKA),I12,IBLTP,NBLOCK,IBLKFO,
     &             I_AM_OUT,N_ELIMINATED_BATCHES)
*.Flush local memory
      CALL GETMEM('KLJ   ','FREE','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','FREE','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','FREE','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXB  ','FREE','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLH1D ','FREE','REAL',KLH1D ,NACOB)
      CALL GETMEM('KLASTR','FREE','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','FREE','INTE',KLBSTR,MXNSTR*NBEL)
      CALL GETMEM('KLRJKA','FREE','REAL',KLRJKA,MAXA)
*
*
C?    STOP ' Jeppe forced me to stop after GASDIA '
      RETURN
      END
