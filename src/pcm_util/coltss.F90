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

subroutine ColTss(IOut,Colour,C1,C2,C3)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: IOut
character(len=*), intent(in) :: Colour
real(kind=wp), intent(out) :: C1, C2, C3

! Assign tesserae colours for GeomView:

if (Colour == 'White') then
  C1 = 1.0_wp
  C2 = 1.0_wp
  C3 = 1.0_wp
else if (Colour == 'Gray') then
  C1 = 0.66_wp
  C2 = 0.66_wp
  C3 = 0.66_wp
else if ((Colour == 'Blue') .or. (Colour == 'Dark Blue')) then
  C1 = 0.0_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if (Colour == 'Light Blue') then
  C1 = 0.0_wp
  C2 = 1.0_wp
  C3 = 1.0_wp
else if (Colour == 'Green') then
  C1 = 0.0_wp
  C2 = 1.0_wp
  C3 = 0.0_wp
else if (Colour == 'Yellow') then
  C1 = 1.0_wp
  C2 = 1.0_wp
  C3 = 0.0_wp
else if (Colour == 'Orange') then
  C1 = 1.0_wp
  C2 = 0.5_wp
  C3 = 0.0_wp
else if (Colour == 'Violet') then
  C1 = 0.6_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if ((Colour == 'Pink') .or. (Colour == 'Light Red')) then
  C1 = 1.0_wp
  C2 = 0.5_wp
  C3 = 1.0_wp
else if (Colour == 'Fuchsia') then
  C1 = 1.0_wp
  C2 = 0.0_wp
  C3 = 1.0_wp
else if ((Colour == 'Red') .or. (Colour == 'Dark Red')) then
  C1 = 1.0_wp
  C2 = 0.0_wp
  C3 = 0.0_wp
else if (Colour == 'Black') then
  C1 = 0.0_wp
  C2 = 0.0_wp
  C3 = 0.0_wp
else
  C1 = 0.0_wp  ! dummy assignement
  C2 = 0.0_wp  ! dummy assignement
  C3 = 0.0_wp  ! dummy assignement
  write(IOut,'(a)') 'Unrecognized colour in ColTss'
  call Abend()
end if

return

end subroutine ColTss
