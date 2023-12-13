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

subroutine basis2run()

use Basis_Info, only: dbsc, iCnttp_Dummy, nCnttp, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iAng, iAtoms, iBasis, icnt, iCnttp, iCo, index_center, iPrim, iShSrt, iyy, jSh, kExp, mdc, nPrim
integer(kind=iwp), allocatable :: IndC(:), primitive_ids(:,:)
real(kind=wp), allocatable :: primitives(:,:)

iAtoms = 0
!***********************************************************************
! Generate list of primitive basis functions
!
! Loop over distinct shell types
nPrim = 0
! Loop over basis sets
do iCnttp=1,nCnttp
  if (iCnttp == iCnttp_Dummy) cycle
  if (dbsc(iCnttp)%iVal == 0) cycle
  mdc = dbsc(iCnttp)%mdci
  iShSrt = dbsc(iCnttp)%iVal
  ! Loop over distinct centers
  do icnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    ! Loop over symmetry-related centers
    do iCo=0,nIrrep/dc(mdc)%nStab-1
      ! Loop over shells associated with this center
      ! Start with s type shells
      jSh = iShSrt
      if (Shells(jSh)%Aux .or. Shells(jSh)%Frag) cycle
      do iAng=0,dbsc(iCnttp)%nVal-1
        nPrim = nPrim+Shells(jSh)%nExp*Shells(jSh)%nBasis
        jSh = jSh+1
      end do
    end do
  end do
end do

call put_iScalar('nPrim',nPrim)

call mma_allocate(IndC,2*S%mCentr,label='IndC')
call mma_allocate(primitive_ids,3,nPrim,label='primitive_ids')
call mma_allocate(primitives,2,nPrim,label='primitives')

! Loop over distinct shell types
iPrim = 0
! Loop over basis sets
do iCnttp=1,nCnttp
  if (iCnttp == iCnttp_Dummy) cycle
  if (dbsc(iCnttp)%iVal == 0) cycle
  mdc = dbsc(iCnttp)%mdci
  iShSrt = dbsc(iCnttp)%iVal
  ! Loop over distinct centers
  do icnt=1,dbsc(iCnttp)%nCntr
    mdc = mdc+1
    ! Loop over symmetry-related centers
    do iCo=0,nIrrep/dc(mdc)%nStab-1
      ! Loop over shells associated with this center
      ! Start with s type shells
      jSh = iShSrt
      if (Shells(jSh)%Aux .or. ShellS(jSh)%Frag) cycle
      ! Get the flat, desymmetrized id of the center
      iyy = Index_Center(mdc,iCo,IndC,iAtoms,S%mCentr)
      do iAng=0,dbsc(iCnttp)%nVal-1
        ! Pointer to the untouched contraction matrix as after input.
        do iBasis=1,Shells(jSh)%nBasis
          do kExp=1,Shells(jSh)%nExp
            iPrim = iPrim+1
            primitive_ids(1,iPrim) = iyy
            primitive_ids(2,iPrim) = iAng
            primitive_ids(3,iPrim) = iBasis
            primitives(1,iPrim) = Shells(jSh)%Exp(kExp)
            primitives(2,iPrim) = Shells(jSh)%Cff_c(kExp,iBasis,2)
          end do
        end do
        jSh = jSh+1
      end do
    end do
  end do

end do

call put_iArray('primitive ids',primitive_ids,3*nPrim)
call put_dArray('primitives',primitives,2*nPrim)

call mma_deallocate(primitive_ids)
call mma_deallocate(primitives)
call mma_deallocate(IndC)

return

end subroutine basis2run
