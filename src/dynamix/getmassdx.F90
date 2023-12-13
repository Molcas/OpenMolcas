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

subroutine GetMassDx(Mass,natom)

use Isotopes, only: Isotope
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: natom
real(kind=wp), intent(out) :: Mass(natom)
integer(kind=iwp) :: matom, i, Iso
character(len=2), allocatable :: atom(:)

call mma_allocate(atom,natom)
call Get_Name_Full(atom)
call Get_nAtoms_All(matom)
call Get_Mass_All(Mass,matom)
do i=matom+1,natom
  atom(i) = adjustl(atom(i))
  Iso = 0
  call Isotope(Iso,atom(i),Mass(i))
end do
call mma_deallocate(atom)

end subroutine GetMassDx
