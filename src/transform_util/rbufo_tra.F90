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

subroutine RBufO_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS)

implicit real*8(a-h,o-z)
#include "SysDef.fh"
#include "stdalloc.fh"
dimension W(*)

call mma_maxDBLE(MEMX)
IADX = (KKTU-1)*IADXS
IST = 1
Length = LBuf
IEnd = LBuf

52 continue
call dDAFILE(LUHLFX,2,W(IST),Length,IADX)
IST = IST+LBuf
IEnd = IEnd+LBuf
if (IEnd > LL) Length = mod(LL,LBuf)
IADX = IADX+(NOTU-1)*IADXS
if (IST <= LL) GO TO 52

IST = 1
return

end subroutine RBufO_tra2
