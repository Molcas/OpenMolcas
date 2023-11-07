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

subroutine CofMss(Coor,nsAtom,cMass)
!***********************************************************************
!   Object: To calculate the molecular mass, the center of mass and    *
!           move the coordinates so origo is the center of mass.       *
!***********************************************************************

use Slapaf_Info, only: dMass, Smmtrc
use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nsAtom
real(kind=wp), intent(in) :: Coor(3,nsAtom)
real(kind=wp), intent(out) :: cMass(3)
integer(kind=iwp) :: i, iCOM, j
real(kind=wp) :: TMass
integer(kind=iwp), external :: iDeg

! Calculate the molecular mass.

TMass = Zero
do I=1,nsAtom
  TMass = TMass+dMass(I)*real(iDeg(Coor(1,i)),kind=wp)
end do
iCOM = -1
if (TMass >= 1.0e99_wp) then
  do i=1,nsAtom
    if (dMass(i) == 1.0e99_wp) then
      iCOM = i
      exit
    end if
  end do
end if

! calculate the center of mass

cMass(:) = Zero
! Loop over the unique centers
do i=1,nsAtom
  do j=1,3
    ! Add contribution
    if (Smmtrc(j,i)) cMass(j) = cMass(j)+dMass(i)*Coor(j,i)*real(iDeg(Coor(1,i)),kind=wp)
  end do
end do

cMass(:) = cMass(:)/TMass
if ((iCOM >= 1) .and. (iCOM <= nsAtom)) cMass(:) = Coor(:,iCom)

#ifdef _DEBUGPRINT_
write(u6,100) cMass(:),TMass
100 format(//,' Center of Mass (Bohr) ',3F10.5,/,' Molecular Mass   (au) ',1F15.5)
#endif

#ifdef _DO_NOT_
! translate the center of mass to origo

do i=1,nsAtom
  Coor(:,i) = Coor(:,i)-cMass(:)
end do
#endif

return

end subroutine CofMss
