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

subroutine Get_Prim_Atom_Tab(nAtoms,nPrim,Work,CENTX,CENTY,CENTZ)

use MPProp_globals, only: iAtPrTab, nAtomPBas
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nPrim
real(kind=wp), intent(in) :: Work(3,nAtoms), CENTX(nPrim*(nPrim+1)/2), CENTY(nPrim*(nPrim+1)/2), CENTZ(nPrim*(nPrim+1)/2)
! EB real(kind=wp), intent(in) :: CENTX(nPrim), CENTY(nPrim), CENTZ(nPrim)
integer(kind=iwp) :: i, j, jj
real(kind=wp), parameter :: Thr = 1.0e-10_wp

!---- Get the primitive basis that belongs to a specific atom

do i=1,nAtoms
  nAtomPBas(i) = 0
  do j=1,nPrim
    jj = j*(j+1)/2
    if ((abs(Work(1,i)-CENTX(jj)) <= Thr) .and. (abs(Work(2,i)-CENTY(jj)) <= Thr) .and. (abs(Work(3,i)-CENTZ(jj)) <= Thr)) then
      nAtomPBas(i) = nAtomPBas(i)+1
      iAtPrTab(nAtomPBas(i),i) = j
      !write(u6,*) 'Atom tab',i,j,CENTX(jj),CENTY(jj),CENTZ(jj),nAtomPBas(i)
    end if
  end do
end do

return

end subroutine Get_Prim_Atom_Tab
