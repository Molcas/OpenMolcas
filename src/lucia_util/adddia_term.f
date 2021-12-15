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
* Copyright (C) 2011, Jeppe Olsen                                      *
*               2011, Giovanni Li Manni                                *
************************************************************************
      SUBROUTINE ADDDIA_TERM(FACTOR,CVEC,SVEC,IASPGP,IBSPGP,IASM,IBSM)
*. Update Sigma vector with diagonal terms for a given block
*     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
*                         + (FACTOR+DIAG(IASPGP,IBSPGP))CVEC(IASPGP,IBSPGP)
*
* Jeppe Olsen and Giovanni Li Manni, September 2011
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "WrkSpc.fh"
#include "mxpdim.fh"
#include "strinp.fh"
#include "orbinp.fh"
#include "gasstr.fh"
#include "strbas.fh"
#include "csm.fh"
#include "stinf.fh"
#include "cecore.fh"
#include "cprnt.fh"
*

*#include <cstate.fh>




*. Input
      DIMENSION CVEC(*)
*. Output
      DIMENSION SVEC(*)
*
      IDUM=0
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA)
*
      IATP   = 1
      IBTP   = 2
      NAEL   = NELEC(IATP)
      NBEL   = NELEC(IBTP)
      NOCTPA = NOCTYP(IATP)
*
      IF(NTEST.GE.10) THEN
        WRITE(6,*) ' ========================='
        WRITE(6,*) '   ADDDIA_TERM for BK is speaking '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' NAEL NBEL =', NAEL,NBEL
        WRITE(6,*) ' IASPGP, IBSPGP = ', IASPGP, IBSPGP
      END IF
*. A bit of scracth
      CALL GETMEM('KLH1D ','ALLO','REAL',KLH1D ,NTOOB)
      CALL GETMEM('KLJ   ','ALLO','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','ALLO','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','ALLO','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXA  ','ALLO','REAL',KLXA  ,NACOB)
      CALL GETMEM('KLXB  ','ALLO','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLSCR ','ALLO','REAL',KLSCR ,2*NACOB)
*. Space for blocks of strings
      CALL GETMEM('KLASTR','ALLO','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','ALLO','INTE',KLBSTR,MXNSTR*NBEL)

      MAXA = IMNMX(IWORK(KNSTSO(IATP)),NSMST*NOCTPA,2)
      CALL GETMEM('KLRJKA','ALLO','REAL',KLRJKA,MAXA)
*. Diagonal of one-body integrals and coulomb and exchange integrals
*. Integrals assumed in place so :
      CALL GT1DIA(WORK(KLH1D))
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS,IREOST)
*. Core energy not included
      ECOREP = 0.0D0
      CALL GTJK(WORK(KLJ),WORK(KLK),NTOOB,WORK(KLSCR2),IREOTS,IREOST)
*
      SHIFT = ECORE_ORIG-ECORE
      FACTORX = FACTOR + SHIFT
*
      CALL ADDDIA_TERMS(NAEL,IWORK(KLASTR),NBEL,IWORK(KLBSTR),
     &             NACOB,CVEC,SVEC,NSMST,WORK(KLH1D),
     &             WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &             WORK(KLK),IWORK(KNSTSO(IATP)),IWORK(KNSTSO(IBTP)),
     &             ECOREP,
     &             IPRDIA,NTOOB,
     &             WORK(KLRJKA),
     &             IASPGP,IASM,IBSPGP,IBSM,FACTORX)
*.Flush local memory
      CALL GETMEM('KLH1D ','FREE','REAL',KLH1D ,NTOOB)
      CALL GETMEM('KLJ   ','FREE','REAL',KLJ   ,NTOOB**2)
      CALL GETMEM('KLK   ','FREE','REAL',KLK   ,NTOOB**2)
      CALL GETMEM('KLSC2 ','FREE','REAL',KLSCR2,2*NTOOB**2)
      CALL GETMEM('KLXA  ','FREE','REAL',KLXA  ,NACOB)
      CALL GETMEM('KLXB  ','FREE','REAL',KLXB  ,NACOB)
      CALL GETMEM('KLSCR ','FREE','REAL',KLSCR ,2*NACOB)
      CALL GETMEM('KLASTR','FREE','INTE',KLASTR,MXNSTR*NAEL)
      CALL GETMEM('KLBSTR','FREE','INTE',KLBSTR,MXNSTR*NBEL)
      CALL GETMEM('KLRJKA','FREE','REAL',KLRJKA,MAXA)

      RETURN
      END
