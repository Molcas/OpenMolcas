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
! Copyright (C) Roland Lindh                                           *
!               2018, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Get_Mass_All
!
!> @brief
!>   Get atomic masses from RUNFILE
!> @author Ignacio Fdez. Galn&aacute;n
!>
!> @details
!> Place atomic masses (in a.u.) into array \p Mass_All(*).
!>
!> @param[out] Mass_All   Array of masses
!> @param[in]  nAtoms_All Number of atoms
!***********************************************************************

subroutine Get_Mass_All(Mass_All,nAtoms_All)

use Symmetry_Info, only: iOper, nIrrep, Symmetry_Info_Get
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms_All
real(kind=wp), intent(out) :: Mass_All(nAtoms_All)
integer(kind=iwp) :: Active = 0, i, iChAtom, iCo, iCoSet(0:7,0:7), iGen(3), iStab(0:7), j, MaxDCR, nAtoms, nAtoms_Allx, nCoSet, &
                     nGen, nStab
real(kind=wp), allocatable :: CU(:,:), Mass(:)
integer(kind=iwp), external :: iChxyz

if (Active == 0) then
  call Symmetry_Info_Get()
  Active = 1
end if
! Obtain symmetry-unique masses
call Get_nAtoms_All(nAtoms_Allx)
if (nAtoms_All /= nAtoms_Allx) then
  write(u6,*) 'Get_Coord_All: nAtoms_All /= nAtoms_Allx'
  write(u6,*) 'nAtoms_All=',nAtoms_All
  write(u6,*) 'nAtoms_Allx=',nAtoms_Allx
  call Abend()
end if
call Get_iScalar('Unique atoms',nAtoms)
call mma_allocate(Mass,nAtoms)
call Get_Mass(Mass,nAtoms)

! Replicate masses
call mma_allocate(CU,3,nAtoms)
call Get_dArray('Unique Coordinates',CU,3*nAtoms)
nGen = 0
if (nIrrep == 2) nGen = 1
if (nIrrep == 4) nGen = 2
if (nIrrep == 8) nGen = 3
if (nGen >= 1) iGen(1) = iOper(1)
if (nGen >= 2) iGen(2) = iOper(2)
if (nGen == 3) iGen(3) = iOper(4)
MaxDCR = 0
j = 0
do i=1,nAtoms
  iChAtom = iChxyz(CU(:,i),iGen,nGen)
  call Stblz(iChAtom,nStab,iStab,MaxDCR,iCoSet)
  nCoSet = nIrrep/nStab
  do iCo=0,nCoSet-1
    j = j+1
    Mass_all(j) = Mass(i)
  end do
end do
call mma_deallocate(CU)
call mma_deallocate(Mass)

return

end subroutine Get_Mass_All
