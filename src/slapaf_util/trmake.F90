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

use Slapaf_Info, only: dMass, Smmtrc
use Symmetry_Info, only: VarR, VarT

implicit real*8(a-h,o-z)
#include "real.fh"
#include "print.fh"
real*8 TRVec(6,3*nAtoms), Coor(3,nAtoms), uMtrx(3*nAtoms), CM(3)
logical SymDsp, CofM

iRout = 131
iPrint = nPrint(iRout)
if (iPrint >= 99) then
  call RecPrt(' In TRMake: Coor',' ',Coor,3,nAtoms)
  write(6,*) ' nDim=',nDim
end if

call dcopy_(6*3*nAtoms,[Zero],0,TRVec,1)
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
      call dcopy_(nAtoms,[One],0,TRVec(nTR,i),18)
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
      j = (iAtom-1)*3+i
      if (CofM) then
        rNorm = rNorm+uMtrx(j)*dMass(iAtom)
        if (Smmtrc(i,iAtom)) then
          CM(i) = CM(i)+uMtrx(j)*Coor(i,iAtom)*dMass(iAtom)
        end if
      else
        rNorm = rNorm+uMtrx(j)
        if (Smmtrc(i,iAtom)) then
          CM(i) = CM(i)+uMtrx(j)*Coor(i,iAtom)
        end if
      end if
    end do
    CM(i) = CM(i)/rNorm
  end do
  !write(6,*) 'TrMake CM=',CM

  do i=1,3
    j = i+1
    if (j > 3) j = j-3
    k = i+2
    if (k > 3) k = k-3
    ! j and k are the index of the plane perpendicular to the axis

    ! Check the rotation has any mirror plane parallel to the axis
    ! of the rotation. If not then the rotation will not break the
    ! symmetry.

    iCmp = 2**(j-1)+2**(k-1)
    if (SymDsp(iCmp)) then
      nTR = nTR+1
      call DYaX(nAtoms,One,Coor(j,1),3,TRVec(nTR,k),18)
      call DaXpY_(nAtoms,-One,CM(j),0,TRVec(nTR,k),18)
      call DYaX(nAtoms,-One,Coor(k,1),3,TRVec(nTR,j),18)
      call DaXpY_(nAtoms,One,CM(k),0,TRVec(nTR,j),18)
    end if
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Normalize vectors

do i=1,nTR
  rii = Zero
  do iAtom=1,3*nAtoms
    rii = rii+uMtrx(iAtom)*TRVec(i,iatom)**2
  end do
  if (rii > 1.d-15) then
    call DScal_(3*nAtoms,One/sqrt(rii),TRVec(i,1),6)
  else
    call dcopy_(3*nAtoms,[Zero],0,TRVec(i,1),6)
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
