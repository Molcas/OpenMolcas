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

subroutine RBufO_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IADXS)

use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: LUHLFX, LL, LBuf, NOTU, KKTU, IADXS
real(kind=wp), intent(_OUT_) :: W(*)
integer(kind=iwp) :: IADX, IEnd, IST, Length, MEMX

call mma_maxDBLE(MEMX)
IADX = (KKTU-1)*IADXS
IST = 1
Length = LBuf
IEnd = LBuf

do
  call dDAFILE(LUHLFX,2,W(IST),Length,IADX)
  IST = IST+LBuf
  IEnd = IEnd+LBuf
  if (IEnd > LL) Length = mod(LL,LBuf)
  IADX = IADX+(NOTU-1)*IADXS
  if (IST > LL) exit
end do

return

end subroutine RBufO_tra2
