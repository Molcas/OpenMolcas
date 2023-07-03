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

implicit real*8(a-h,o-z)
logical Torsion_Check
real*8 Ref(3,4), Coor(3,4)
integer iTabAtoms(2,0:nMax,mAtoms), iCase(4,12), iSet(4), jSet(4)
data iCase/1,2,3,4,1,2,4,3,2,1,3,4,2,1,4,3,1,3,2,4,1,3,4,2,3,1,2,4,3,1,4,2,1,4,2,3,1,4,3,2,4,1,2,3,4,1,3,2/

Torsion_Check = .false.
iSet(1) = iAtom
iSet(2) = jAtom
iSet(3) = kAtom
iSet(4) = lAtom
!write(6,*) 'Here we go!'
!write(6,*) 'iSet=',iSet

! Check all twelve possible torsions to see which one is the one
! which should at least be included.

FC = 0.0d0
Test = 0.0d0
do ii=1,12
  jSet(1) = iSet(iCase(1,ii))
  jSet(2) = iSet(iCase(2,ii))
  jSet(3) = iSet(iCase(3,ii))
  jSet(4) = iSet(iCase(4,ii))
  !write(6,*) 'ii,jSet=',ii,jSet
  call dcopy_(3,Ref(1,iCase(1,ii)),1,Coor(1,1),1)
  call dcopy_(3,Ref(1,iCase(2,ii)),1,Coor(1,2),1)
  call dcopy_(3,Ref(1,iCase(3,ii)),1,Coor(1,3),1)
  call dcopy_(3,Ref(1,iCase(4,ii)),1,Coor(1,4),1)
  !call RecPrt('Coor',' ',Coor,3,4)
  Test = FC_Torsion(jSet,Coor,iTabAtoms,nMax,mAtoms)
  if (Test > FC) then
    if (ii > 1) return
    FC = Test
  end if
end do
Torsion_Check = .true.
!write(6,*) 'Torsion_Check=',Torsion_Check

return

end function Torsion_Check
