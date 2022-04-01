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

subroutine Cooout(Head,Cordst,nPart,nCent)

use Definitions, only: wp, iwp, u6

implicit none
character(len=200), intent(in) :: Head
integer(kind=iwp), intent(in) :: nPart, nCent
real(kind=wp), intent(in) :: Cordst(3,nPart*nCent)
integer(kind=iwp) :: i, j, kaunter

write(u6,*)
write(u6,*)
write(u6,'(A)') Head
kaunter = 0
do i=1,nPart
  write(u6,*) 'Molecule ',i
  do j=1,nCent
    kaunter = kaunter+1
    write(u6,*) Cordst(:,kaunter)
  end do
end do

return

end subroutine Cooout
