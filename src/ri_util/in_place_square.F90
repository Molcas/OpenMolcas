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

subroutine In_place_Square(Buff,nBuff)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nBuff
real(kind=wp), intent(inout) :: Buff(nBuff,nBuff)
integer(kind=iwp) :: i, j

!call RecPrt('Buff',' ',Buff,nBuff,nBuff)
do j=1,nBuff
  do i=1,j-1
    Buff(j,i) = Buff(i,j)
  end do
end do
!call RecPrt('Buff',' ',Buff,nBuff,nBuff)
!write(u6,'(10F10.3)') (Buff(i,i),i=1,nBuff)

return

end subroutine In_place_Square
