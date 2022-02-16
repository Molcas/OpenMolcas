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

use mrci_global, only: CSPCK, FOCK, IFIRST, INTSY, IRC, ISAB, ISMAX, KBUFF1, NBITM1, NBITM2, NBITM3, NBTRI, NCHN1, NCHN2, NCHN3, &
                       NVIRT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NHDIAG, NINTGR, NVT
integer(kind=iwp), allocatable :: Inds(:,:)
real(kind=wp), allocatable :: ACBDS(:), ACBDT(:), BACBD(:), BIAC1(:), BICA1(:), BUFBI(:), Bufs(:,:), FIIJJ(:), FIJIJ(:), HDIAG(:)

! ----------------------------------------------------------------------
call mma_allocate(Bufs,NBITM1,NCHN1,label='Bufs')
call mma_allocate(Inds,NBITM1+2,NCHN1,label='Inds')
call mma_allocate(BUFBI,KBUFF1,label='BUFBI')
call mma_allocate(BIAC1,ISMAX,label='BIAC1')
call mma_allocate(BICA1,ISMAX,label='BICA1')
Bufs(:,:) = Zero
Inds(:,:) = 0
call SORTA(Bufs,Inds,ISAB,BUFBI,BIAC1,BICA1,NINTGR)
call mma_deallocate(Bufs)
call mma_deallocate(Inds)
call mma_deallocate(BUFBI)
call mma_deallocate(BIAC1)
call mma_deallocate(BICA1)
! ----------------------------------------------------------------------

if (IFIRST == 0) then
  call mma_allocate(Bufs,NBITM2,NCHN2,label='Bufs')
  call mma_allocate(Inds,NBITM2+2,NCHN2,label='Bufs')
  call mma_allocate(BACBD,KBUFF1,label='BACBD')
  call mma_allocate(ACBDT,ISMAX,label='ACBDT')
  call mma_allocate(ACBDS,ISMAX,label='ACBDS')
  Bufs(:,:) = Zero
  Inds(:,:) = 0
  call SORTB(Bufs,Inds,ACBDS,ACBDT,ISAB,BACBD)
  call mma_deallocate(Bufs)
  call mma_deallocate(Inds)
  call mma_deallocate(BACBD)
  call mma_deallocate(ACBDT)
  call mma_deallocate(ACBDS)
end if
! ------------------- SORT --------------------------------------------
call mma_allocate(Bufs,NBITM3,NCHN3,label='Bufs')
call mma_allocate(Inds,NBITM3+2,NCHN3,label='Inds')
call mma_allocate(FIIJJ,NBTRI,label='FIIJJ')
call mma_allocate(FIJIJ,NBTRI,label='FIJIJ')
Bufs(:,:) = Zero
Inds(:,:) = 0
call SORT_MRCI(Bufs,Inds,FOCK,FIIJJ,FIJIJ)
call mma_deallocate(Bufs)
call mma_deallocate(Inds)
! ----------------------------------------------------------------------
NVT = (NVIRT*(NVIRT+1))/2
NHDIAG = max(NVT,IRC(1))
call mma_allocate(HDIAG,NHDIAG,label='HDIAG')
call IIJJ(CSPCK,INTSY,HDIAG,FOCK,FIIJJ,FIJIJ)
call IJIJ(INTSY,HDIAG,FIJIJ)
call mma_deallocate(FIIJJ)
call mma_deallocate(FIJIJ)
call mma_deallocate(HDIAG)

return

end subroutine DIAGCT
