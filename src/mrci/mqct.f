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
      SUBROUTINE MQCT(AREF,EREF,CI,SGM,ICI)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "SysDef.fh"

#include "mrci.fh"
#include "WrkSpc.fh"
      DIMENSION AREF(NREF,NREF),EREF(NREF)
      DIMENSION ICI(MBUF),SGM(NCONF),CI(NCONF)
*      LCBUF =LMQ
*      LSBUF =LCBUF +MXVEC*MBUF
*      LDBUF =LSBUF +MXVEC*MBUF
*      LCSECT=LDBUF +MBUF
*      LRSECT=LCSECT+NSECT*MXVEC
*      LXI1  =LRSECT+NSECT*MXVEC
*      LXI2  =LXI1  +NSECT*NRROOT
*      LCNEW =LXI2  +NSECT*NRROOT
*      LSCR  =LCNEW +NSECT*NRROOT
      CALL GETMEM('CBUF','ALLO','REAL',LCBUF,MXVEC*MBUF)
      CALL GETMEM('SBUF','ALLO','REAL',LSBUF,MXVEC*MBUF)
      CALL GETMEM('DBUF','ALLO','REAL',LDBUF,MBUF)
      CALL GETMEM('CSECT','ALLO','REAL',LCSECT,NSECT*MXVEC)
      CALL GETMEM('RSECT','ALLO','REAL',LRSECT,NSECT*MXVEC)
      CALL GETMEM('XI1','ALLO','REAL',LXI1,NSECT*NRROOT)
      CALL GETMEM('XI2','ALLO','REAL',LXI2,NSECT*NRROOT)
      CALL GETMEM('CNEW','ALLO','REAL',LCNEW,NSECT*NRROOT)
      CALL GETMEM('SCR','ALLO','REAL',LSCR,MXVEC*MBUF)
      CALL DIAGRO(CI,SGM,Work(LCBUF),
     *            Work(LSBUF),Work(LDBUF),
     *            AREF,EREF,
     *            Work(LCSECT),Work(LRSECT),
     *            Work(LXI1),Work(LXI2),Work(LCNEW),
     *            Work(LSCR),ICI)
      CALL GETMEM('CBUF','FREE','REAL',LCBUF,MXVEC*MBUF)
      CALL GETMEM('SBUF','FREE','REAL',LSBUF,MXVEC*MBUF)
      CALL GETMEM('DBUF','FREE','REAL',LDBUF,MBUF)
      CALL GETMEM('CSECT','FREE','REAL',LCSECT,NSECT*MXVEC)
      CALL GETMEM('RSECT','FREE','REAL',LRSECT,NSECT*MXVEC)
      CALL GETMEM('XI1','FREE','REAL',LXI1,NSECT*NRROOT)
      CALL GETMEM('XI2','FREE','REAL',LXI2,NSECT*NRROOT)
      CALL GETMEM('CNEW','FREE','REAL',LCNEW,NSECT*NRROOT)
      CALL GETMEM('SCR','FREE','REAL',LSCR,MXVEC*MBUF)
      RETURN
      END
