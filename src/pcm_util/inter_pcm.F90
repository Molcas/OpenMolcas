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

subroutine Inter_PCM(XE,YE,ZE,RE,P1,P2,P3,P4,I,IPRINT)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: XE, YE, ZE, RE, P1(3), P2(3), P3(3)
real(kind=wp), intent(out) :: P4(3)
integer(kind=iwp), intent(in) :: I, IPRINT
integer(kind=iwp) :: M
real(kind=wp) :: ALPHA, DELTA, DIFF, DIFF2, DNORM, R, R2
real(kind=wp), parameter :: TOL = 1.0e-12_wp

! Trova il punto P4, sull'arco P1-P2 sotteso dal centro P3, che
! si trova sulla superficie della sfera XE,YE,ZE,RE
! P4 e' definito come combinazioe lineare di P1 e P2, con
! il parametro ALPHA ottimizzato per tentativi.

R2 = (P1(1)-P3(1))**2+(P1(2)-P3(2))**2+(P1(3)-P3(3))**2
R = sqrt(R2)
ALPHA = Half
DELTA = Zero
M = 1
do
  if (M > 100) then
    if (IPRINT > 0) write(u6,'(/,10X," INTER: too many iterations")')
    return
  end if
  ALPHA = ALPHA+DELTA
  DNORM = Zero
  P4(:) = P1(:)+ALPHA*(P2(:)-P1(:))-P3(:)
  DNORM = DNORM+P4(1)*P4(1)+P4(2)*P4(2)+P4(3)*P4(3)
  DNORM = sqrt(DNORM)
  P4(:) = P4(:)*R/DNORM+P3(:)
  DIFF2 = (P4(1)-XE)**2+(P4(2)-YE)**2+(P4(3)-ZE)**2
  DIFF = sqrt(DIFF2)-RE
  if (abs(DIFF) < TOL) exit
  if (I == 0) then
    if (DIFF > Zero) DELTA = One/(Two**(M+1))
    if (DIFF < Zero) DELTA = -One/(Two**(M+1))
    M = M+1
  else
    if (DIFF > Zero) DELTA = -One/(Two**(M+1))
    if (DIFF < Zero) DELTA = One/(Two**(M+1))
    M = M+1
  end if
end do

return

end subroutine Inter_PCM
