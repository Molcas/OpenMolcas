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

subroutine RBuf_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUHLFX, LL, LBuf, NOTU, KKTU, IADXS
real(kind=wp), intent(_OUT_) :: W(*)
integer(kind=iwp), intent(out) :: IST
!integer(kind=iwp) :: BLKSZ, MEMX, NBLCK

!call mma_maxDBLE(MEMX)
!BLKSZ = (NOTU-1)*IADXS
!NBLCK = MEMX/BLKSZ

!if ((NBLCK > 1) .and. (LBuf < 256)) then
!  call RBufF_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IADXS,MEMX)
!else
call RBufO_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IADXS)
!end if
IST = 1

return

end subroutine RBuf_tra2
