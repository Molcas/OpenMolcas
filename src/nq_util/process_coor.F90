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

subroutine Process_Coor(R,Coor,nAtoms,nSym,iOper)

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: R(3)
real(kind=wp), intent(inout) :: Coor(3,*)
integer(kind=iwp), intent(inout) :: nAtoms
integer(kind=iwp), intent(in) :: nSym, iOper(0:nSym-1)
integer(kind=iwp) :: iAtom, iRef, iSym
real(kind=wp) :: Q(3)

!                                                                      *
!***********************************************************************
!                                                                      *
!call RecPrt('Coor(Enter)',' ',Coor,3,nAtoms)
!call RecPrt('R',' ',R,3,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Identify if this is a new center.

do iAtom=1,nAtoms
  if ((R(1) == Coor(1,iAtom)) .and. (R(2) == Coor(2,iAtom)) .and. (R(3) == Coor(3,iAtom))) return
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Add this atom to the list

nAtoms = nAtoms+1
Coor(:,nAtoms) = R
iRef = nAtoms
!call RecPrt('Coor(updated)',' ',Coor,3,nAtoms)
!write(u6,*) 'nSym=',nSym
!write(u6,*) 'iOper=',iOper
!                                                                      *
!***********************************************************************
!                                                                      *
! Add symmetry degenerate atoms to the list

outer: do iSym=1,nSym-1
  !write(u6,*) 'iOper(iSym)=',iOper(iSym)
  Q(:) = R
  if (btest(iOper(iSym),0)) Q(1) = -Q(1)
  if (btest(iOper(iSym),1)) Q(2) = -Q(2)
  if (btest(iOper(iSym),2)) Q(3) = -Q(3)
  !call RecPrt('Q',' ',Q,3,1)
  do iAtom=iRef,nAtoms
    if ((Q(1) == Coor(1,iAtom)) .and. (Q(2) == Coor(2,iAtom)) .and. (Q(3) == Coor(3,iAtom))) cycle outer
  end do
  nAtoms = nAtoms+1
  Coor(:,nAtoms) = Q
end do outer
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine Process_Coor
