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

subroutine BornMayerBK(iQ_Atoms,BoMaH,BoMaO)
! With the Brdarski-Karlstrom scheme, construct the Born-Mayer parameters.

use qmstat_global, only: CharDi, CharDiQ, iPrint, QuaDi, QuaDiQ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iQ_Atoms
real(kind=wp), intent(out) :: BoMaH(iQ_Atoms), BoMaO(iQ_Atoms)
integer(kind=iwp) :: i, j
real(kind=wp) :: rBdi(2), rdi2
real(kind=wp), allocatable :: rBdiQ(:)
real(kind=wp), parameter :: cjhr = 0.1734_wp ! What is this number?

! The solvent part.

do i=1,2
  rdi2 = Zero
  do j=1,3
    rdi2 = rdi2+QuaDi(j,i)
  end do
  rBdi(i) = sqrt(rdi2/CharDi(i))
end do

! The solute part.

call mma_allocate(rBdiQ,iQ_Atoms,label='rBdiQ')
do i=1,iQ_Atoms
  rdi2 = Zero
  do j=1,3
    rdi2 = rdi2+QuaDiQ(j,i)
  end do
  rBdiQ(i) = sqrt(rdi2/CharDiQ(i))
end do

! Put together.

BoMaH(:) = One/(cjhr*(RBdiQ(:)+rbdi(1)))
BoMaO(:) = One/(cjhr*(RBdiQ(:)+rbdi(2)))
if (iPrint >= 8) then
  write(6,*) '   Born-Mayer parameters.'
  do i=1,iQ_Atoms
    write(6,'(A,i2,A,2(f12.4))') '    Atom ',i,' (H/O):',BoMaH(i),BoMaO(i)
  end do
end if
call mma_deallocate(rBdiQ)

return

end subroutine BornMayerBK
