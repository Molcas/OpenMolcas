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

subroutine dRdR(NSI,NESFJ,DR,RSolv,Sphere,NewSph)

use Constants, only: Two, Three, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSI, NESFJ, NewSph(2)
real(kind=wp), intent(out) :: DR
real(kind=wp), intent(in) :: RSolv, Sphere(4,*)
integer(kind=iwp) :: NSK
real(kind=wp) :: D, D2, RI, RJ, RK, RS

! Trova la derivata del raggio della sfera NSI rispetto al raggio
! della sfera NSJ.
!
! La sfera NSI (che appartiene alle sfere "aggiunte" da PEDRA)
! dipende dalle due sfere "precedenti" NESFJ e NSK
! Se NESFJ o NSK sono negativi, la sfera aggiunta e' di tipo C
! e la generatrice "principale" corrisponde al label negativo
! (cfr. JCC 11, 1047 (1990))

if ((NEWSPH(1) < 0) .or. (NEWSPH(2) < 0)) then
  NSK = NEWSPH(1)
  if (abs(NSK) == NESFJ) NSK = NEWSPH(2)
  if (NSK > 0) then
    RS = RSOLV
    RJ = Sphere(4,NESFJ)+RS
    RK = Sphere(4,NSK)+RS
    RI = Sphere(4,NSI)+RS
    D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
    D = sqrt(D2)
    DR = (Two*D*RJ+Two*D*Sphere(4,NESFJ)-Two*RJ*Sphere(4,NESFJ)+D*D-RJ*RJ-RK*RK)/(Two*D*RI)
  else
    RS = RSOLV
    RJ = Sphere(4,NESFJ)+RS
    RI = Sphere(4,NSI)+RS
    D2 = (Sphere(1,NESFJ)-Sphere(1,abs(NSK)))**2+(Sphere(2,NESFJ)-Sphere(2,abs(NSK)))**2+(Sphere(3,NESFJ)-Sphere(3,abs(NSK)))**2
    D = sqrt(D2)
    DR = (Sphere(4,abs(NSK))*RJ)/(D*RI)
  end if
else
  NSK = NEWSPH(1)
  if (NSK == NESFJ) NSK = NEWSPH(2)
  RS = RSOLV
  RJ = Sphere(4,NESFJ)+RS
  RK = Sphere(4,NSK)+RS
  RI = Sphere(4,NSI)+RS
  D2 = (Sphere(1,NESFJ)-Sphere(1,NSK))**2+(Sphere(2,NESFJ)-Sphere(2,NSK))**2+(Sphere(3,NESFJ)-Sphere(3,NSK))**2
  D = sqrt(D2)
  DR = (-Three*RJ*RJ+RK*RK+Two*RJ*RK+Three*D*RJ-D*RK)/(Four*D*RI)
end if

return

end subroutine dRdR
