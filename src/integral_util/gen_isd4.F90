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

subroutine Gen_iSD4(iS,jS,kS,lS,iSD,nSD,iSD4)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iS, jS, kS, lS, nSD, iSD(0:nSD,1024)
integer(kind=iwp), intent(out) :: iSD4(0:nSD,4)
integer(kind=iwp) :: i, iSkal, j, jQuad(4)

jQuad(1) = iS
jQuad(2) = jS
jQuad(3) = kS
jQuad(4) = lS
do i=1,4
  iSkal = jQuad(i)
  do j=0,nSD
    iSD4(j,i) = iSD(j,iSkal)
  end do
end do

return

end subroutine Gen_iSD4
