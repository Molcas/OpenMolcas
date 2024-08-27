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

subroutine Stblz(iChxyz,nStab,jStab,MaxDCR,iCoSet)

use Symmetry_Info, only: iOper, nIrrep

implicit none
integer iChxyz, nStab, MaxDCR
integer jStab(0:7), iCoSet(0:7,0:7)
integer iStab, i, j, nMax, iTest, ielem, iTmp, iOpMn, ip

!                                                                      *
!***********************************************************************
!                                                                      *
! Find the Stabilizers of this center
!                                                                      *
!***********************************************************************
!                                                                      *
iStab = 0
do i=0,nIrrep-1
  if (iand(iChxyz,iOper(i)) == 0) then
    iStab = iStab+1
    jStab(iStab-1) = iOper(i)
  end if
end do
nStab = iStab
MaxDCR = max(MaxDCR,iStab)
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate all possible (left) CoSet
!                                                                      *
!***********************************************************************
!                                                                      *
do i=0,nIrrep-1
  do j=0,iStab-1
    iCoSet(i,j) = ieor(iOper(i),jStab(j))
  end do
end do
if (iStab == 1) Go To 39 ! skip if all are unique

! Order the Coset so we will have the unique ones first

nMax = 1
if (nMax == nIrrep/iStab) Go To 39 ! skip if there is only one
iTest = iStab-1 ! Test on the last element
do j=1,nIrrep-1 !
  ! Check uniqueness
  do i=0,nMax-1
    do ielem=0,iStab-1
      if (iCoSet(i,iTest) == iCoSet(j,ielem)) Go To 35
    end do
  end do
  ! Move unique CoSet to nMax+1
  nMax = nMax+1
  do ielem=0,iStab-1
    iTmp = iCoSet(nMax-1,ielem)
    iCoSet(nMax-1,ielem) = iCoSet(j,ielem)
    iCoSet(j,ielem) = iTmp
  end do
  if (nMax == nIrrep/iStab) Go To 39
35 continue
end do
39 continue
!                                                                      *
!***********************************************************************
!                                                                      *
do i=0,nIrrep/iStab-1
  iOpMn = iCoSet(i,0)
  do j=1,iStab-1
    iOpMn = iand(iOpMn,iCoset(i,j))
  end do
  ip = 0
  do j=0,iStab-1
    if (iOpMn == iCoSet(i,j)) ip = j
  end do
  iTmp = iCoSet(i,0)
  iCoSet(i,0) = iCoSet(i,ip)
  iCoSet(i,ip) = iTmp

end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Stblz
