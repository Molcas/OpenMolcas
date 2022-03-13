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

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iS, nSt, iDm
real(kind=wp) :: Dens(*), StVec(iDm,*)
integer(kind=iwp) :: i, ii, j, jj, kaunt

kaunt = 0
do i=1,nSt
  do j=1,i
    kaunt = kaunt+1
    Dens(kaunt) = Zero
  end do
end do
kaunt = 0
do ii=1,nSt
  do jj=1,ii
    kaunt = kaunt+1
    if (ii == jj) then
      Dens(kaunt) = StVec(ii,iS)*StVec(jj,iS)
    else
      Dens(kaunt) = Two*StVec(ii,iS)*StVec(jj,iS)
    end if
  end do
end do

return

end subroutine DensiSt
