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

subroutine BPut(EVec,nDim,BMx,nX,Smmtrc,nQQ,Degen)

implicit real*8(a-h,o-z)
#include "real.fh"
real*8 EVec(nDim,nQQ), BMx(nX,nQQ), Degen(nX)
logical Smmtrc(nX)

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
iDim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    iDim = iDim+1
    do iQQ=1,nQQ
      BMx(iX,iQQ) = EVec(iDim,iQQ)/sqrt(Degen(iX))
    end do
  else
    do iQQ=1,nDim
      BMx(iX,iQQ) = Zero
    end do
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('BPut: BMx',' ',BMx,nX,nQQ)
#endif

return

end subroutine BPut
