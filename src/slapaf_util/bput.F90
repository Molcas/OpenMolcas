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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nDim, nX, nQQ
real(kind=wp) :: EVec(nDim,nQQ), BMx(nX,nQQ), Degen(nX)
logical(kind=iwp) :: Smmtrc(nX)
integer(kind=iwp) :: i_Dim, iQQ, iX

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
i_Dim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    i_Dim = i_Dim+1
    do iQQ=1,nQQ
      BMx(iX,iQQ) = EVec(i_Dim,iQQ)/sqrt(Degen(iX))
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
