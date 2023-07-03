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

subroutine ShfANM(nInter,nIter,rInt,Shift)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 rInt(nInter,nIter), Shift(nInter,nIter)

if (nIter == 1) return

#ifdef _DEBUGPRINT_
call RecPrt(' ShfANM: rInt',' ',rInt,nInter,nIter)
#endif

! Shifts: dq = q   -q
!           n   n+1  n

do Iter=1,nIter-1
  call dcopy_(nInter,rInt(1,Iter+1),1,Shift(1,Iter),1)
  call DaXpY_(nInter,-One,rInt(1,Iter),1,Shift(1,Iter),1)
end do

#ifdef _DEBUGPRINT_
call RecPrt(' In ShfANM: New Shifts',' ',Shift,nInter,nIter-1)
#endif

return

end subroutine ShfANM
