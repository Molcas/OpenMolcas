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

subroutine SqToTri_Q(SqMat,TriMat,iDi)

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iDi
real(kind=wp), intent(in) :: SqMat(iDi,iDi)
real(kind=wp), intent(out) :: TriMat(nTri_Elem(iDi))
integer(kind=iwp) :: i, j, kaunter

kaunter = 0
do i=1,iDi
  do j=1,i
    kaunter = kaunter+1
    TriMat(kaunter) = SqMat(i,j)
  end do
end do

return

end subroutine SqToTri_Q
