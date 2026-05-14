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
integer(kind=iwp), intent(in) :: nDim, nX, nQQ
real(kind=wp), intent(in) :: EVec(nDim,nQQ), Degen(nX)
real(kind=wp), intent(out) :: BMx(nX,nQQ)
logical(kind=iwp) :: Smmtrc(nX)
integer(kind=iwp) :: i_Dim, iX

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
i_Dim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    i_Dim = i_Dim+1
    BMx(iX,:) = EVec(i_Dim,:)/sqrt(Degen(iX))
  else
    BMx(iX,:) = Zero
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('BPut: BMx',' ',BMx,nX,nQQ)
#endif

return

end subroutine BPut
