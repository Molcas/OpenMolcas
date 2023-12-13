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

subroutine PrePro(nAtom,Coor)

use Slapaf_Info, only: Atom, Grd, iInt, iRow, iter, lNmHss, lOld, lOld_Implicit, mRowH, mTROld, nDimBC, nFix, nSup, Redundant
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtom
real(kind=wp), intent(in) :: Coor(3,nAtom)
integer(kind=iwp) :: mTR, nQQ, nRowH
logical(kind=iwp) :: CofM
real(kind=wp), allocatable :: TR(:)

CofM = (Iter == 1) .and. lNmHss
call mma_allocate(TR,18*nAtom,Label='TR')
TR(:) = Zero
call TRPGen(nDimBC,nAtom,Coor,mTR,CofM,TR)
call mma_deallocate(TR)
if (lNmHss) then
  if (Iter == 1) mTROld = mTR
  if ((iter <= 2*(nDimBC-mTROld)+1) .and. (iter /= 1)) mTR = mTROld
else
  mTROld = mTR
end if

! Operate according to two modes
! iRow > 0 : user supplied internal coordinates
! iRow <= 0 : Cartesian or Internal Coordinates

nRowH = 0
if (iRow > 0) then

  ! Find the number of active and frozen internal coordinates.

  call Rd_UDIC(iInt,nFix,nRowH)
  nQQ = iInt+nFix
  if (nRowH > 0) then
    call mma_allocate(mRowH,nRowH,Label='mRowH')
    call Rd_UDIC_RowH(nQQ,nRowH,mRowH)
  end if
  if (nQQ > nDimBC-mTR) Redundant = .true.

else

  nFix = 0
end if

! Initiate the force constant matrix in internal coordinate
! basis, excluding rotation and translation.
! Write to runfile only on the first iteration and that there
! was not an already defined Hessian.

if (iter == 1) call IntFcm(lOld_Implicit)
if ((.not. lOld) .and. lOld_Implicit) lOld = .true.

! Symmetrize forces

if (allocated(nSup)) then
  call SupSym(Grd,nAtom,Coor,size(nSup),nSup,Atom)
  call mma_deallocate(Atom)
  call mma_deallocate(nSup)
end if

return

end subroutine PrePro
