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
      SUBROUTINE DIAGCT()
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
#include "WrkSpc.fh"
      ILIM=4
      IF(IFIRST.NE.0)ILIM=2
* ----------------------------------------------------------------------
      NINDS=(NBITM1+2)*NCHN1
      NBUFS=NBITM1*NCHN1
      CALL GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Allo','Inte',LINDS,NINDS)
      NBUFBI=KBUFF1
      CALL GETMEM('BUFBI','Allo','Real',LBUFBI,NBUFBI)
      CALL GETMEM('BIAC1','Allo','Real',LBIAC1,ISMAX)
      CALL GETMEM('BICA1','Allo','Real',LBICA1,ISMAX)
      CALL DCOPY_(NBUFS,[0.0D0],0,Work(LBUFS),1)
      CALL ICOPY(NINDS,[0],0,IWork(LINDS),1)
      CALL SORTA (Work(LBUFS),IWork(LINDS),IWork(LISAB),
     *            Work(LBUFBI),Work(LBIAC1),Work(LBICA1),NINTGR)
      CALL GETMEM('BIAC1','Free','Real',LBIAC1,ISMAX)
      CALL GETMEM('BICA1','Free','Real',LBICA1,ISMAX)
      CALL GETMEM('BUFBI','Free','Real',LBUFBI,NBUFBI)
      CALL GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Free','Inte',LINDS,NINDS)
* ----------------------------------------------------------------------

      IF(IFIRST.EQ.0) THEN
      NINDS=(NBITM2+2)*NCHN2
      NBUFS=NBITM2*NCHN2
      CALL GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Allo','Inte',LINDS,NINDS)
      CALL GETMEM('BACBD','Allo','Real',LBACBD,KBUFF1)
      CALL GETMEM('ACBDT','Allo','Real',LACBDT,ISMAX)
      CALL GETMEM('ACBDS','Allo','Real',LACBDS,ISMAX)
      CALL DCOPY_(NBUFS,[0.0D0],0,Work(LBUFS),1)
      CALL ICOPY(NINDS,[0],0,IWork(LINDS),1)
      CALL SORTB (Work(LBUFS),IWork(LINDS),
     *            Work(LACBDS),Work(LACBDT),IWork(LISAB),
     *            Work(LBACBD),NINTGR)
      CALL GETMEM('BACBD','Free','Real',LBACBD,KBUFF1)
      CALL GETMEM('ACBDT','Free','Real',LACBDT,ISMAX)
      CALL GETMEM('ACBDS','Free','Real',LACBDS,ISMAX)
      CALL GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Free','Inte',LINDS,NINDS)
      END IF
* ------------------- SORT --------------------------------------------
      NINDS=(NBITM3+2)*NCHN3
      NBUFS=NBITM3*NCHN3
      CALL GETMEM('Bufs','Allo','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Allo','Inte',LINDS,NINDS)
      CALL GETMEM('FIIJJ','Allo','Real',LIIJJ,NBTRI)
      CALL GETMEM('FIJIJ','Allo','Real',LIJIJ,NBTRI)
      CALL DCOPY_(NBUFS,[0.0D0],0,Work(LBUFS),1)
      CALL ICOPY(NINDS,[0],0,IWork(LINDS),1)
      CALL SORT_MRCI (Work(LBUFS),IWork(LINDS),Work(LFOCK),Work(LIIJJ),
     *           Work(LIJIJ),NINTGR)
      CALL GETMEM('Bufs','Free','Real',LBUFS,NBUFS)
      CALL GETMEM('Inds','Free','Inte',LINDS,NINDS)
* ----------------------------------------------------------------------
      NVT=(NVIRT*(NVIRT+1))/2
      NHDIAG=MAX(NVT,IRC(1))
      CALL GETMEM('HDIAG','Allo','Real',LHDIAG,NHDIAG)
      CALL IIJJ (IWork(LCSPCK),IWork(LINTSY),Work(LHDIAG),
     *           Work(LFOCK),Work(LIIJJ),Work(LIJIJ))
      CALL GETMEM('FIIJJ','Free','Real',LIIJJ,NBTRI)
      CALL IJIJ (IWork(LINTSY),Work(LHDIAG),Work(LFOCK),Work(LIJIJ))
      CALL GETMEM('HDIAG','Free','Real',LHDIAG,NHDIAG)
      CALL GETMEM('FIJIJ','Free','Real',LIJIJ,NBTRI)
      RETURN
      END
