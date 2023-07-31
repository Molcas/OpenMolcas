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

function Torsion_Check(iAtom,jAtom,kAtom,lAtom,Ref,iTabAtoms,nMax,mAtoms)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp) :: Torsion_Check
integer(kind=iwp), intent(in) :: iAtom, jAtom, kAtom, lAtom, nMax, mAtoms, iTabAtoms(2,0:nMax,mAtoms)
real(kind=wp), intent(in) :: Ref(3,4)
integer(kind=iwp) :: ii, iSet(4), jSet(4)
real(kind=wp) :: Coor(3,4), FC, Test
integer(kind=iwp), parameter :: iCase(4,12) = reshape([1,2,3,4, &
                                                       1,2,4,3, &
                                                       2,1,3,4, &
                                                       2,1,4,3, &
                                                       1,3,2,4, &
                                                       1,3,4,2, &
                                                       3,1,2,4, &
                                                       3,1,4,2, &
                                                       1,4,2,3, &
                                                       1,4,3,2, &
                                                       4,1,2,3, &
                                                       4,1,3,2],[4,12])
real(kind=wp), external :: FC_Torsion

Torsion_Check = .false.
iSet(1) = iAtom
iSet(2) = jAtom
iSet(3) = kAtom
iSet(4) = lAtom
!write(u6,*) 'Here we go!'
!write(u6,*) 'iSet=',iSet

! Check all twelve possible torsions to see which one is the one
! which should at least be included.

FC = Zero
Test = Zero
do ii=1,12
  jSet(1) = iSet(iCase(1,ii))
  jSet(2) = iSet(iCase(2,ii))
  jSet(3) = iSet(iCase(3,ii))
  jSet(4) = iSet(iCase(4,ii))
  !write(u6,*) 'ii,jSet=',ii,jSet
  Coor(:,1) = Ref(:,iCase(1,ii))
  Coor(:,2) = Ref(:,iCase(2,ii))
  Coor(:,3) = Ref(:,iCase(3,ii))
  Coor(:,4) = Ref(:,iCase(4,ii))
  !call RecPrt('Coor',' ',Coor,3,4)
  Test = FC_Torsion(jSet,Coor,iTabAtoms,nMax,mAtoms)
  if (Test > FC) then
    if (ii > 1) return
    FC = Test
  end if
end do
Torsion_Check = .true.
!write(u6,*) 'Torsion_Check=',Torsion_Check

return

end function Torsion_Check
