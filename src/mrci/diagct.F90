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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "mrci.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: NBUFBI, NHDIAG, NINTGR, NVT
integer(kind=iwp), allocatable :: Inds(:,:)
real(kind=wp), allocatable :: Bufs(:,:)

! ----------------------------------------------------------------------
call mma_allocate(Bufs,NBITM1,NCHN1,label='Bufs')
call mma_allocate(Inds,NBITM1+2,NCHN1,label='Inds')
NBUFBI = KBUFF1
call GETMEM('BUFBI','Allo','Real',LBUFBI,NBUFBI)
call GETMEM('BIAC1','Allo','Real',LBIAC1,ISMAX)
call GETMEM('BICA1','Allo','Real',LBICA1,ISMAX)
Bufs(:,:) = Zero
Inds(:,:) = 0
call SORTA(Bufs,Inds,IWork(LISAB),Work(LBUFBI),Work(LBIAC1),Work(LBICA1),NINTGR)
call GETMEM('BIAC1','Free','Real',LBIAC1,ISMAX)
call GETMEM('BICA1','Free','Real',LBICA1,ISMAX)
call GETMEM('BUFBI','Free','Real',LBUFBI,NBUFBI)
call mma_deallocate(Bufs)
call mma_deallocate(Inds)
! ----------------------------------------------------------------------

if (IFIRST == 0) then
  call mma_allocate(Bufs,NBITM2,NCHN2,label='Bufs')
  call mma_allocate(Inds,NBITM2+2,NCHN2,label='Bufs')
  call GETMEM('BACBD','Allo','Real',LBACBD,KBUFF1)
  call GETMEM('ACBDT','Allo','Real',LACBDT,ISMAX)
  call GETMEM('ACBDS','Allo','Real',LACBDS,ISMAX)
  Bufs(:,:) = Zero
  Inds(:,:) = 0
  call SORTB(Bufs,Inds,Work(LACBDS),Work(LACBDT),IWork(LISAB),Work(LBACBD),NINTGR)
  call GETMEM('BACBD','Free','Real',LBACBD,KBUFF1)
  call GETMEM('ACBDT','Free','Real',LACBDT,ISMAX)
  call GETMEM('ACBDS','Free','Real',LACBDS,ISMAX)
  call mma_deallocate(Bufs)
  call mma_deallocate(Inds)
end if
! ------------------- SORT --------------------------------------------
call mma_allocate(Bufs,NBITM3,NCHN3,label='Bufs')
call mma_allocate(Inds,NBITM3+2,NCHN3,label='Inds')
call GETMEM('FIIJJ','Allo','Real',LIIJJ,NBTRI)
call GETMEM('FIJIJ','Allo','Real',LIJIJ,NBTRI)
Bufs(:,:) = Zero
Inds(:,:) = 0
call SORT_MRCI(Bufs,Inds,Work(LFOCK),Work(LIIJJ),Work(LIJIJ),NINTGR)
call mma_deallocate(Bufs)
call mma_deallocate(Inds)
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
