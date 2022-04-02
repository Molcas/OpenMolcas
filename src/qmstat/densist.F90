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

! The RASSI-density matrix subroutine.
subroutine DensiSt(Dens,StVec,iS,nSt,iDm)
! iS    - Which state that is occupied.
! Dens  - The density
! StVec - The coefficients for how the new states are expressed with the old.

use Index_Functions, only: nTri_Elem
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iS, nSt, iDm
real(kind=wp), intent(out) :: Dens(nTri_Elem(nSt))
real(kind=wp), intent(in) :: StVec(iDm,*)
integer(kind=iwp) :: i, j, kaunt

kaunt = 0
do i=1,nSt
  do j=1,i-1
    kaunt = kaunt+1
    Dens(kaunt) = Two*StVec(i,iS)*StVec(j,iS)
  end do
  kaunt = kaunt+1
  Dens(kaunt) = StVec(i,iS)*StVec(i,iS)
end do

return

end subroutine DensiSt
