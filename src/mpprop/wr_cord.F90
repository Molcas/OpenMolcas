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

subroutine Wr_Cord(nAtoms)

use MPProp_globals, only: Cor, iAtomType, iAtomPar, Labe, Qnuc
use Constants, only: Angstrom
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms
integer(kind=iwp) :: i, iStdOut

iStdOut = u6

write(iStdOut,*)
write(iStdOut,'(10X,A)') ' ************************************************ '
write(iStdOut,'(10X,A)') ' **** Cartesian Coordinates / Bohr, Angstrom **** '
write(iStdOut,'(10X,A)') ' ************************************************ '
write(iStdOut,*)
write(iStdOut,'(A11,A7,A)') 'Center','Label', &
                            '              X              Y              Z                     X              Y              Z'
do i=1,nAtoms
  !Jose write(iStdOut,'(6X,I3,5X,A4,3(5X,F10.6),7X,3(5X,F10.6))')
  write(iStdOut,'(6X,I3,5X,A6,3(5X,F10.6),7X,3(5X,F10.6))') i,Labe(i),Cor(:,i,i),Cor(:,i,i)*Angstrom
end do
write(iStdOut,*)
write(iStdOut,*)
write(iStdOut,'(A11,A7,A13,A25,A20)') 'Center','Label','Atomtype','Effective core charge','Atomic Parameter'
write(iStdOut,*)
do i=1,nAtoms
  !Jose write(iStdOut,'(6X,I3,5X,A4,10X,I3,19X,F6.1,10X,I10)') i,Labe(i),
  write(iStdOut,'(6X,I3,5X,A6,10X,I3,19X,F6.1,10X,I10)') i,Labe(i),iAtomType(i),Qnuc(i),iAtomPar(i)
end do
write(iStdOut,*)

return

end subroutine Wr_Cord
