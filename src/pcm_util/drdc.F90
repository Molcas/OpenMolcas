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

subroutine dRdC(NSI,ICOORD,NESFJ,DR,RSolv,Sphere,NewSph)

use Constants, only: Two, Four, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSI, ICOORD, NESFJ, NewSph(2)
real(kind=wp), intent(out) :: DR
real(kind=wp), intent(in) :: RSolv, Sphere(4,*)
integer(kind=iwp) :: K, NSK
real(kind=wp) :: A, B, COORDJ(3), COORDK(3), D, D2, DIFF, FAC, RI, RJ, RK, RS

! Trova la derivata del raggio della sfera NSI rispetto alla
! coordinata ICOORD (1=X, 2=Y, 3=Z) della sfera NSJ, che interseca
! NSI.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e K
!
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1) < 0) .or. (NEWSPH(2) < 0)) then
  NSK = NEWSPH(1)
  if (abs(NSK) == NESFJ) NSK = NEWSPH(2)
  COORDJ(:) = Sphere(1:3,NESFJ)
  COORDK(:) = Sphere(1:3,abs(NSK))
  RI = Sphere(4,NSI)+RSOLV
  RJ = Sphere(4,NESFJ)+RSOLV
  RK = Sphere(4,abs(NSK))+RSOLV
  DIFF = COORDJ(ICOORD)-COORDK(ICOORD)
  D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
  D = sqrt(D2)
  FAC = Sphere(4,NESFJ)*(RJ*RJ-D*D-RK*RK)
  if (NSK < 0) FAC = Sphere(4,abs(NSK))*(RK*RK-D*D-RJ*RJ)
  DR = DIFF*FAC/(Two*D**3*RI)
else
  K = NEWSPH(1)
  if (K == NESFJ) K = NEWSPH(2)
  COORDJ(:) = Sphere(1:3,NESFJ)
  COORDK(:) = Sphere(1:3,K)
  D2 = (COORDJ(1)-COORDK(1))**2+(COORDJ(2)-COORDK(2))**2+(COORDJ(3)-COORDK(3))**2
  D = sqrt(D2)
  B = Half*(D+Sphere(4,NESFJ)-Sphere(4,K))
  RS = RSOLV
  A = ((Sphere(4,NESFJ)+RS)**2+D2-(Sphere(4,K)+RS)**2)/D
  DR = (Two*A*B-Two*B*D-A*D)*(COORDJ(ICOORD)-COORDK(ICOORD))/(Four*D2*(Sphere(4,NSI)+RS))
end if

return

end subroutine dRdC
