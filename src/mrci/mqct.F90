!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine MQCT(AREF,EREF,CI,SGM,ICI)

use Definitions, only: wp, iwp

implicit none
#include "mrci.fh"
real(kind=wp) :: AREF(NREF,NREF), EREF(NREF), CI(NCONF), SGM(NCONF)
integer(kind=iwp) :: ICI(MBUF)
#include "WrkSpc.fh"
integer(kind=iwp) :: LCBUF, LCNEW, LCSECT, LDBUF, LRSECT, LSBUF, LSCR, LXI1, LXI2

!LCBUF = LMQ
!LSBUF = LCBUF+MXVEC*MBUF
!LDBUF = LSBUF+MXVEC*MBUF
!LCSECT = LDBUF+MBUF
!LRSECT = LCSECT+NSECT*MXVEC
!LXI1 = LRSECT+NSECT*MXVEC
!LXI2 = LXI1+NSECT*NRROOT
!LCNEW = LXI2+NSECT*NRROOT
!LSCR = LCNEW+NSECT*NRROOT
call GETMEM('CBUF','ALLO','REAL',LCBUF,MXVEC*MBUF)
call GETMEM('SBUF','ALLO','REAL',LSBUF,MXVEC*MBUF)
call GETMEM('DBUF','ALLO','REAL',LDBUF,MBUF)
call GETMEM('CSECT','ALLO','REAL',LCSECT,NSECT*MXVEC)
call GETMEM('RSECT','ALLO','REAL',LRSECT,NSECT*MXVEC)
call GETMEM('XI1','ALLO','REAL',LXI1,NSECT*NRROOT)
call GETMEM('XI2','ALLO','REAL',LXI2,NSECT*NRROOT)
call GETMEM('CNEW','ALLO','REAL',LCNEW,NSECT*NRROOT)
call GETMEM('SCR','ALLO','REAL',LSCR,MXVEC*MBUF)
call DIAGRO(CI,SGM,Work(LCBUF),Work(LSBUF),Work(LDBUF),AREF,EREF,Work(LCSECT),Work(LRSECT),Work(LXI1),Work(LXI2),Work(LCNEW), &
            Work(LSCR),ICI)
call GETMEM('CBUF','FREE','REAL',LCBUF,MXVEC*MBUF)
call GETMEM('SBUF','FREE','REAL',LSBUF,MXVEC*MBUF)
call GETMEM('DBUF','FREE','REAL',LDBUF,MBUF)
call GETMEM('CSECT','FREE','REAL',LCSECT,NSECT*MXVEC)
call GETMEM('RSECT','FREE','REAL',LRSECT,NSECT*MXVEC)
call GETMEM('XI1','FREE','REAL',LXI1,NSECT*NRROOT)
call GETMEM('XI2','FREE','REAL',LXI2,NSECT*NRROOT)
call GETMEM('CNEW','FREE','REAL',LCNEW,NSECT*NRROOT)
call GETMEM('SCR','FREE','REAL',LSCR,MXVEC*MBUF)

return

end subroutine MQCT
