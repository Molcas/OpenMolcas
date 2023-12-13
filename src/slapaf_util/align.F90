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
! Copyright (C) 2013, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Align
!
!> @brief
!>   Align two structures.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Align a molecular structure with the reference, using the stored weights
!> (masses by default). This is sometimes needed to ensure that an optimal
!> structure is found when there are constraints expressed in weighted space
!> (e.g. `sphere` or `transverse`).
!>
!> @param[in,out] Coord Cartesian coordinates to align
!> @param[in]     Ref   Cartesian coordinates of the reference structure
!> @param[in]     nAtom Number of symmetry-unique atoms
!***********************************************************************

subroutine Align(Coord,Ref,nAtom)

use Symmetry_Info, only: iOper, nIrrep, VarR, VarT
use Slapaf_Info, only: Weights
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), intent(inout) :: Coord(3*nAtom)
real(kind=wp), intent(in) :: Ref(3*nAtom)
integer(kind=iwp) :: i, iAdr, iAt, iChxyz, iIrrep, mAtom, nStb
real(kind=wp) :: RMS, RMSMax
integer(kind=iwp), allocatable :: iStab(:)
real(kind=wp), allocatable :: Coor_All(:,:), Ref_All(:,:)

! Do nothing if the energy is not rot. and trans. invariant
if (VarR .or. VarT) return

call mma_allocate(Coor_All,3,nAtom*8,Label='Coor_All')
call Expand_Coor(Coord,nAtom,Coor_All,mAtom)
call mma_allocate(Ref_All,3,nAtom*8,Label='Ref_All')
call Expand_Coor(Ref,nAtom,Ref_All,mAtom)

!call RecPrt('Coord before align',' ',Coor_All,3,mAtom)

call Superpose_w(Coor_All,Ref_All,Weights,mAtom,RMS,RMSMax)

! Get the stabilizers for each atom (to keep the symmetry)
! (code copied from init_slapaf)
call mma_allocate(iStab,nAtom,Label='iStab')
do iAt=1,nAtom
  iAdr = (iAt-1)*3+1
  iChxyz = 0
  do i=0,2
    if (Ref(iAdr+i) /= Zero) then
      do iIrrep=0,nIrrep-1
        if (btest(iOper(iIrrep),i)) iChxyz = ibset(iChxyz,i)
      end do
    end if
  end do
  nStb = 0
  do iIrrep=0,nIrrep-1
    if ((nStb <= 1) .and. (iand(iChxyz,iOper(iIrrep)) == 0)) then
      iStab(iAt) = iOper(iIrrep)
      nStb = nStb+1
    end if
  end do
end do

call Fix_Symmetry(Coor_All,nAtom,iStab)
call mma_deallocate(iStab)

Coord(:) = reshape(Coor_All(:,1:nAtom),[3*nAtom])

!call RecPrt('Coord after align',' ',Coor_All,3,mAtom)

call mma_deallocate(Coor_All)
call mma_deallocate(Ref_All)

end subroutine Align
