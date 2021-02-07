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

use Isotopes

implicit none
#include "stdalloc.fh"
real*8, intent(inout) :: Mass(*)
integer, intent(in) :: natom
integer :: matom, i, Iso
character, allocatable :: atom(:)*2

call mma_allocate(atom,natom)
call Get_Name_Full(atom)
call Get_nAtoms_All(matom)
call Get_Mass_All(Mass,matom)
do i=1,natom
  if (i > matom) then
    call LeftAd(atom(i))
    Iso = 0
    call Isotope(Iso,atom(i),Mass(i))
  end if
end do
call mma_deallocate(atom)

end subroutine GetMassDx
