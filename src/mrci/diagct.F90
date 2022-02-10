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

subroutine DIAGCT()

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"

! ----------------------------------------------------------------------
NINDS = (NBITM1+2)*NCHN1
NBUFS = NBITM1*NCHN1
call GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
call GETMEM('Inds','Allo','Inte',LINDS,NINDS)
NBUFBI = KBUFF1
call GETMEM('BUFBI','Allo','Real',LBUFBI,NBUFBI)
call GETMEM('BIAC1','Allo','Real',LBIAC1,ISMAX)
call GETMEM('BICA1','Allo','Real',LBICA1,ISMAX)
call DCOPY_(NBUFS,[0.0d0],0,Work(LBUFS),1)
call ICOPY(NINDS,[0],0,IWork(LINDS),1)
call SORTA(Work(LBUFS),IWork(LINDS),IWork(LISAB),Work(LBUFBI),Work(LBIAC1),Work(LBICA1),NINTGR)
call GETMEM('BIAC1','Free','Real',LBIAC1,ISMAX)
call GETMEM('BICA1','Free','Real',LBICA1,ISMAX)
call GETMEM('BUFBI','Free','Real',LBUFBI,NBUFBI)
call GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
call GETMEM('Inds','Free','Inte',LINDS,NINDS)
! ----------------------------------------------------------------------

if (IFIRST == 0) then
  NINDS = (NBITM2+2)*NCHN2
  NBUFS = NBITM2*NCHN2
  call GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
  call GETMEM('Inds','Allo','Inte',LINDS,NINDS)
  call GETMEM('BACBD','Allo','Real',LBACBD,KBUFF1)
  call GETMEM('ACBDT','Allo','Real',LACBDT,ISMAX)
  call GETMEM('ACBDS','Allo','Real',LACBDS,ISMAX)
  call DCOPY_(NBUFS,[0.0d0],0,Work(LBUFS),1)
  call ICOPY(NINDS,[0],0,IWork(LINDS),1)
  call SORTB(Work(LBUFS),IWork(LINDS),Work(LACBDS),Work(LACBDT),IWork(LISAB),Work(LBACBD),NINTGR)
  call GETMEM('BACBD','Free','Real',LBACBD,KBUFF1)
  call GETMEM('ACBDT','Free','Real',LACBDT,ISMAX)
  call GETMEM('ACBDS','Free','Real',LACBDS,ISMAX)
  call GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
  call GETMEM('Inds','Free','Inte',LINDS,NINDS)
end if
! ------------------- SORT --------------------------------------------
NINDS = (NBITM3+2)*NCHN3
NBUFS = NBITM3*NCHN3
call GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
call GETMEM('Inds','Allo','Inte',LINDS,NINDS)
call GETMEM('FIIJJ','Allo','Real',LIIJJ,NBTRI)
call GETMEM('FIJIJ','Allo','Real',LIJIJ,NBTRI)
call DCOPY_(NBUFS,[0.0d0],0,Work(LBUFS),1)
call ICOPY(NINDS,[0],0,IWork(LINDS),1)
call SORT_MRCI(Work(LBUFS),IWork(LINDS),Work(LFOCK),Work(LIIJJ),Work(LIJIJ),NINTGR)
call GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
call GETMEM('Inds','Free','Inte',LINDS,NINDS)
! ----------------------------------------------------------------------
NVT = (NVIRT*(NVIRT+1))/2
NHDIAG = max(NVT,IRC(1))
call GETMEM('HDIAG','Allo','Real',LHDIAG,NHDIAG)
call IIJJ(IWork(LCSPCK),IWork(LINTSY),Work(LHDIAG),Work(LFOCK),Work(LIIJJ),Work(LIJIJ))
call GETMEM('FIIJJ','Free','Real',LIIJJ,NBTRI)
call IJIJ(IWork(LINTSY),Work(LHDIAG),Work(LFOCK),Work(LIJIJ))
call GETMEM('HDIAG','Free','Real',LHDIAG,NHDIAG)
call GETMEM('FIJIJ','Free','Real',LIJIJ,NBTRI)

return

end subroutine DIAGCT
