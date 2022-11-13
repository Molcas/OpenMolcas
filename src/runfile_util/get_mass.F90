!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2018, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Get_Mass
!
!> @brief
!>   Get (symmetry-unique) atomic masses from RUNFILE
!> @author Ignacio Fdez. Galn&aacute;n
!>
!> @details
!> Place atomic masses (in a.u.) into array \p Mass_All(*).
!>
!> @param[out] Mass   Array of masses
!> @param[in]  nAtoms Number of atoms
!***********************************************************************

subroutine Get_Mass(Mass,nAtoms)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms
real(kind=wp), intent(out) :: Mass(nAtoms)
integer(kind=iwp) :: i, mAtoms, nCent
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: AtoB(:)
real(kind=wp), allocatable :: CentMass(:)

call Get_iScalar('Unique atoms',mAtoms)
if (mAtoms /= nAtoms) then
  write(u6,*) 'Get_Mass: mAtoms /= nAtoms'
  write(u6,*) 'mAtoms=',mAtoms
  write(u6,*) 'nAtoms=',nAtoms
  call Abend()
end if
call mma_allocate(AtoB,nAtoms)
call Get_iArray('Atom -> Basis',AtoB,nAtoms)
call Qpg_dArray('Isotopes',Found,nCent)
if (.not. Found) then
  write(u6,*) 'Get_Mass: Isotopes array not found'
  call Abend()
end if
call mma_allocate(CentMass,nCent)
call Get_dArray('Isotopes',CentMass,nCent)
do i=1,nAtoms
  Mass(i) = CentMass(AtoB(i))
end do
call mma_deallocate(CentMass)
call mma_deallocate(AtoB)

return

end subroutine Get_Mass
