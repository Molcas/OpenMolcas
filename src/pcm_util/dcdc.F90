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

subroutine dCdC(JJ,ICOORD,NESFJ,DC,Sphere,NewSph)

use Constants, only: One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: JJ, ICOORD, NESFJ, NewSph(2)
real(kind=wp), intent(out) :: DC
real(kind=wp), intent(in) :: Sphere(4,*)
integer(kind=iwp) :: K, NSK
real(kind=wp) :: COORDJ(3), COORDK(3), D, D2

! Trova la derivata della coordinata JJ del centro della sfera
! NewSph rispetto alla coordinata ICOORD di NSJ, che interseca NewSph.
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
  COORDJ(:) = Sphere(1:3,NESFJ)
  COORDK(:) = Sphere(1:3,abs(NSK))
  D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
  D = sqrt(D2)
  if (NSK > 0) then
    DC = Sphere(4,NESFJ)*(COORDJ(JJ)-COORDK(JJ))*(COORDJ(ICOORD)-COORDK(ICOORD))/D**3
    if (ICOORD == JJ) DC = DC+One-Sphere(4,NESFJ)/D
  else
    DC = -Sphere(4,abs(NSK))*(COORDK(JJ)-COORDJ(JJ))*(COORDK(ICOORD)-COORDJ(ICOORD))/D**3
    if (ICOORD == JJ) DC = DC+Sphere(4,abs(NSK))/D
  end if
else
  K = NEWSPH(1)
  if (K == NESFJ) K = NEWSPH(2)
  COORDJ(:) = Sphere(1:3,NESFJ)
  COORDK(:) = Sphere(1:3,K)
  D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
  D = sqrt(D2)
  DC = (Sphere(4,NESFJ)-Sphere(4,K))*(COORDJ(ICOORD)-COORDK(ICOORD))*(COORDJ(JJ)-COORDK(JJ))/(Two*D**3)
  if (JJ == ICOORD) DC = DC+Half-(Sphere(4,NESFJ)-Sphere(4,K))/(Two*D)
end if

return

end subroutine dCdC
