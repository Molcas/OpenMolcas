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

subroutine TRComp(TRVec,nTR,nX,Smmtrc)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTR, nX
real(kind=wp), intent(inout) :: TRVec(nTR,nX)
logical(kind=iwp), intent(in) :: Smmtrc(nX)
integer(kind=iwp) :: i_Dim, iX

if (nTR == 0) return
i_Dim = 0
do iX=1,nX
  if (Smmtrc(iX)) then
    i_Dim = i_Dim+1
    if (i_Dim /= iX) TRVec(:,i_Dim) = TRVec(:,iX)
  end if
end do

return

end subroutine TRComp
