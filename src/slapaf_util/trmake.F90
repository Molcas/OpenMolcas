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

subroutine TRMake(TRVec,Coor,nAtoms,nTR,uMtrx,nDim,CofM)

use Symmetry_Info, only: VarR, VarT
use Slapaf_Info, only: dMass, Smmtrc
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nDim
real(kind=wp), intent(out) :: TRVec(6,3,nAtoms)
real(kind=wp), intent(in) :: Coor(3,nAtoms), uMtrx(3,nAtoms)
integer(kind=iwp), intent(out) :: nTR
logical(kind=iwp), intent(in) :: CofM
#include "print.fh"
integer(kind=iwp) :: i, iAtom, iCmp, iPrint, iRout, j, k
real(kind=wp) :: CM(3), rii, rNorm
logical(kind=iwp) :: SymDsp

iRout = 131
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In TRMake: Coor',' ',Coor,3,nAtoms)
  write(u6,*) ' nDim=',nDim
end if

TRVec(:,:,:) = Zero
nTR = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Translation

if (.not. VarT) then
  do i=1,3
    iCmp = 2**(i-1)
    if (SymDsp(iCmp)) then
      nTR = nTR+1
      TRVec(nTR,i,:) = One

    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Rotation around some center

! Loop over axis

if (.not. VarR) then
  do i=1,3
    CM(i) = Zero
    rNorm = Zero
    do iAtom=1,nAtoms
      if (CofM) then
        rNorm = rNorm+uMtrx(i,iAtom)*dMass(iAtom)
        if (Smmtrc(i,iAtom)) CM(i) = CM(i)+uMtrx(i,iAtom)*Coor(i,iAtom)*dMass(iAtom)
      else
        rNorm = rNorm+uMtrx(i,iAtom)
        if (Smmtrc(i,iAtom)) CM(i) = CM(i)+uMtrx(i,iAtom)*Coor(i,iAtom)
      end if
    end do
    CM(i) = CM(i)/rNorm
  end do
  !write(u6,*) 'TrMake CM=',CM

  do i=1,3
    j = mod(i,3)+1
    k = mod(i+1,3)+1
    ! j and k are the index of the plane perpendicular to the axis

    ! Check the rotation has any mirror plane parallel to the axis
    ! of the rotation. If not then the rotation will not break the
    ! symmetry.

    iCmp = ibset(ibset(0,j-1),k-1)
    if (SymDsp(iCmp)) then
      nTR = nTR+1
      TRVec(nTR,k,:) = Coor(j,:)-CM(j)
      TRVec(nTR,j,:) = CM(k)-Coor(k,:)
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Normalize vectors

do i=1,nTR
  rii = Zero
  do iAtom=1,nAtoms
    do j=1,3
      rii = rii+uMtrx(j,iAtom)*TRVec(i,j,iatom)**2
    end do
  end do
  if (rii > 1.0e-15_wp) then
    TRVec(i,:,:) = TRVec(i,:,:)/sqrt(rii)
  else
    TRVec(i,:,:) = Zero
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 99) call RecPrt(' In TRMake: TRVec',' ',TRVec,6,3*nAtoms)
call TROrder(TRVec,nTR,3*nAtoms)
if (iPrint >= 99) call RecPrt(' In TRMake: TRVec',' ',TRVec,nTR,3*nAtoms)
call TRComp(TRVec,nTR,3*nAtoms,SmmTrc)

if (iPrint >= 99) call RecPrt(' In TRMake: TRVec',' ',TRVec,nTR,nDim)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine TRMake
