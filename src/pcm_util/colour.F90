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

subroutine Colour(NesfP,NAt,AtmC,IAt,Coor_Sph,N,C1,C2,C3)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NesfP, NAt, IAt(*), N
real(kind=wp), intent(in) :: AtmC(3,*), Coor_Sph(4,*)
real(kind=wp), intent(out) :: C1, C2, C3
integer(kind=iwp) :: I, J
real(kind=wp) :: Diff
character(len=20) :: Col
real(kind=wp), parameter :: Delta = 1.0e-3_wp

! Assign tesserae colours for GeomView:
! Carbon: green, nitrogen: blue, oxygen: red, hydrogen: light blue,
! others: violet, added spheres: gray.

Col = ' '
if (N > NESFP) then
  Col = 'Gray'
  call ColTss(u6,Col,C1,C2,C3)
  return
end if
I = 0
do J=1,NAt
  I = I+1
  Diff = sqrt((AtmC(1,J)-Coor_Sph(1,N))**2+(AtmC(2,J)-Coor_Sph(2,N))**2+(AtmC(3,J)-Coor_Sph(3,N))**2)
  if (Diff < Delta) then
    if (IAt(I) == 6) then
      Col = 'Green'
    elseif (IAt(I) == 7) then
      Col = 'Blue'
    elseif (IAt(I) == 8) then
      Col = 'Red'
    elseif (IAt(I) == 1) then
      Col = 'Light Blue'
    else
      Col = 'Fuchsia'
    end if
  end if
end do
call ColTss(u6,Col,C1,C2,C3)

return

end subroutine Colour
