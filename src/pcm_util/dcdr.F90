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

subroutine dCdR(JJ,NESFJ,DC,Sphere,NewSph)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: JJ, NESFJ, NewSph(2)
real(kind=wp), intent(out) :: DC
real(kind=wp), intent(in) :: Sphere(4,*)
integer(kind=iwp) :: NSK
real(kind=wp) :: COORDJ(3), COORDK(3), D, D2

! Trova la derivata della coordinata JJ del centro della sfera
! NewSph rispetto al raggio dellla sfera NSJ.
!
! La sfera NewSph (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e NSK
!
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1) < 0) .or. (NEWSPH(2) < 0)) then
  NSK = NEWSPH(1)
  if (abs(NSK) == NESFJ) NSK = NEWSPH(2)
  DC = Zero
  if (NSK > 0) then
    COORDJ(:) = Sphere(1:3,NESFJ)
    COORDK(:) = Sphere(1:3,NSK)
    D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
    D = sqrt(D2)
    DC = -(COORDJ(JJ)-COORDK(JJ))/D
  end if
else
  NSK = NEWSPH(1)
  if (NSK == NESFJ) NSK = NEWSPH(2)
  COORDJ(:) = Sphere(1:3,NESFJ)
  COORDK(:) = Sphere(1:3,NSK)
  D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
  D = sqrt(D2)
  DC = -(COORDJ(JJ)-COORDK(JJ))/(Two*D)
end if

return

end subroutine dCdR
