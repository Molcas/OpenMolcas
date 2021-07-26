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

subroutine RBufF_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IADXS,MEMX)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUHLFX, LL, LBuf, NOTU, KKTU, IADXS, MEMX
real(kind=wp), intent(_OUT_) :: W(*)
integer(kind=iwp) :: BLKSZ, BPASS, I, IADX, IST, J, NBLCK, NPASS, NRST
real(kind=wp), allocatable :: BUF(:,:)

BLKSZ = (NOTU-1)*IADXS+LBuf
NBLCK = MEMX/BLKSZ

call mma_allocate(BUF,BLKSZ,NBLCK,LABEL='BUF')

BPASS = (LL+LBuf-1)/Lbuf
NPASS = (BPASS+NBLCK-1)/NBLCK

!write(u6,*) 'LL=',LL
!write(u6,*) 'LBUF=',LBUF
!write(u6,*) 'MEMX=',MEMX
!write(u6,*) 'BLKSZ=',BLKSZ
!write(u6,*) 'NBLCK=',NBLCK
!write(u6,*) 'BPASS=',BPASS
!write(u6,*) 'NPASS=',NPASS

IADX = (KKTU-1)*IADXS
IST = 1

do I=1,NPASS-1
  call dDAFILE(LUHLFX,2,BUF,size(BUF),IADX)

  do J=1,NBLCK
    call dcopy_(LBuf,BUF(:,J),1,W(IST),1)
    IST = IST+LBuf
  end do

  BPASS = BPASS-NBLCK
end do

NRST = (BPASS-1)*BLKSZ+mod(LL,LBuf)
call dDAFILE(LUHLFX,2,BUF,NRST,IADX)

do J=1,BPASS-1
  call dcopy_(LBuf,BUF(:,J),1,W(IST),1)
  IST = IST+LBuf
end do

NRST = mod(LL,LBuf)

call dcopy_(NRST,BUF(:,BPASS),1,W(IST),1)

call mma_deallocate(BUF)

return

end subroutine RBufF_tra2
