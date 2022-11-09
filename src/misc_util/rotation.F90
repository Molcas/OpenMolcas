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

subroutine Rotation(TotalM,TRotA,TRotB,TRotC,nsRot,nFAtoms,lSlapaf)

use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half, Angstrom, auTocm, auTokJ, auToHz, kBoltzmann, uToau
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(out) :: TotalM, TRotA, TRotB, TRotC
integer(kind=iwp), intent(inout) :: nsRot
integer(kind=iwp), intent(out) :: nFAtoms
logical(kind=iwp), intent(in) :: lSlapaf
#include "Molcas.fh"
integer(kind=iwp) :: i, iAtom, j, nrot
real(kind=wp) :: CM(3), dEV, dSum, dVec(3), dX, dY, dZ, EVal(nTri_Elem(3)), Inrt(3,3), RotE(3), Vec(3,3)
real(kind=wp), allocatable :: CCoor(:,:), FCoor(:,:), Mass(:), SOCoor(:,:)
character(len=LenIn), allocatable :: FAtLbl(:)
real(kind=wp), parameter :: RT = 1.0e3_wp*Half*auTokJ/kBoltzmann/uToau

!CM              : Center of masses
!FCoor           : Full Coords
!CCoor           : Mass-centered Coords
!SOCoor          : Symmetry-Oriented Coords
!Inrt, RotE, Vec : Inertia components
!Mass            : Masses
!FAtLbl          : Atomic labels

! Get Atomic Full Labels, Coordinates & Mass - FAtLbl, FCoor, Mass

call Get_nAtoms_All(nFAtoms)
call mma_allocate(FCoor,3,nFAtoms,label='FCoor')
call mma_allocate(Mass,nFAtoms,label='Mass')
call mma_allocate(FAtLbl,nFAtoms,label='FAtLbl')
call GetFullCoord(FCoor,Mass,FAtLbl,nFAtoms,lSlapaf)

! Define the Center of Masses - CM()

CM(1) = Zero
CM(2) = Zero
CM(3) = Zero
TotalM = Zero
do i=1,nFAtoms
  TotalM = TotalM+Mass(i)
  CM(:) = CM+Mass(i)*FCoor(:,i)
end do
CM(:) = CM/TotalM

! Shift coordinates in the center of masses - CCoord()

call mma_allocate(CCoor,3,nFAtoms,label='CCoor')
do i=1,nFAtoms
  CCoor(:,i) = FCoor(:,i)-CM
end do
call mma_deallocate(FCoor)

! Compute the Inertia-matrix - Inrt(3,3)

Inrt(:,:) = Zero
do i=1,nFAtoms
  dX = CCoor(1,i)
  dY = CCoor(2,i)
  dZ = CCoor(3,i)
  Inrt(1,1) = Inrt(1,1)+Mass(i)*(dY*dY+dZ*dZ) ! YY ZZ
  Inrt(2,1) = Inrt(2,1)-Mass(i)*dX*dY         ! XY
  Inrt(2,2) = Inrt(2,2)+Mass(i)*(dX*dX+dZ*dZ) ! XX ZZ
  Inrt(3,1) = Inrt(3,1)-Mass(i)*dX*dZ         ! XZ
  Inrt(3,2) = Inrt(3,2)-Mass(i)*dY*dZ         ! YZ
  Inrt(3,3) = Inrt(3,3)+Mass(i)*(dX*dX+dY*dY) ! XX YY
end do
Inrt(1,2) = Inrt(2,1)
Inrt(1,3) = Inrt(3,1)
Inrt(2,3) = Inrt(3,2)

! and diagonalize it - Inrt(3,3)

do i=1,3
  do j=1,i
    EVal(iTri(i,j)) = Inrt(i,j)
  end do
end do
call unitmat(Vec,3)
call Jacob(EVal,Vec,3,3)
call Jacord(EVal,Vec,3,3)
do i=1,3
  RotE(i) = EVal(nTri_Elem(i))
end do

! Sort the principal axis such that z' is the one with the lowest eigenvalue.

do i=1,2
  do j=i+1,3
    if (RotE(i) < RotE(j)) then
      dEV = RotE(i)
      dVec(:) = Vec(:,i)
      RotE(i) = RotE(j)
      Vec(:,i) = Vec(:,j)
      RotE(j) = dEV
      Vec(:,j) = dVec
    end if
  end do
end do

! Rotate coords to Symmetry-Oriented

call mma_allocate(SOCoor,3,nFAtoms,label='SOCoor')
do iAtom=1,nFAtoms
  do i=1,3
    dSum = Zero
    do j=1,3
      dSum = dSum+CCoor(j,iAtom)*Vec(j,i)
    end do
    SOCoor(i,iAtom) = dSum
  end do
end do
call mma_deallocate(CCoor)

! Rotational Symmetry factor - nsRot

if (nsRot == 0) nsRot = 1
if (nFAtoms == 2) then
  if (Mass(1) == Mass(2)) nsRot = 2
end if

TRotA = RT/(RotE(3)+1.0e-99_wp)
TRotB = RT/(RotE(2)+1.0e-99_wp)
TRotC = RT/(RotE(1)+1.0e-99_wp)

! Check if linear molecule

nrot = 3
if (TRotA > 1.0e99_wp) nrot = nrot-1
if (TRotB > 1.0e99_wp) nrot = nrot-1
if (TRotC > 1.0e99_wp) nrot = nrot-1

! Print results

write(u6,'(A)') ' Mass-centered Coordinates (Angstrom):'
write(u6,'(1X,A)') '********************************************************'
write(u6,'(1X,A)') 'Label        X           Y           Z          Mass  '
write(u6,'(1X,A)') '--------------------------------------------------------'
do i=1,nFAtoms
  write(u6,'(1X,A,1X,3F12.6,1x,F12.5)') FAtLbl(i),Angstrom*SOCoor(:,i),Mass(i)
end do
write(u6,'(1X,A)') '--------------------------------------------------------'
write(u6,'(A,F12.6)') ' Molecular mass:',TotalM
write(u6,'(A,3F10.4)') ' Rotational Constants (cm-1):',(auTocm*Half/(uToau*RotE(i)),i=1,nrot)
write(u6,'(A,3F10.4)') ' Rotational Constants (GHz) :',(1.0e-9_wp*auToHz*Half/(uToau*RotE(i)),i=1,nrot)
write(u6,'(A,3F10.4)') ' Rotational temperatures (K):',(RT/RotE(i),i=1,nrot)
write(u6,'(A,I2)') ' Rotational Symmetry factor: ',nsRot
call mma_deallocate(Mass)
call mma_deallocate(FAtLbl)
call mma_deallocate(SOCoor)

return

end subroutine Rotation
