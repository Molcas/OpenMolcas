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

subroutine GauBon(XE,YE,ZE,RE,IntSph,NV,NS,PTS,CCC,PP,AREA,IPRINT)

use PCM_Arrays, only: MxVert
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IntSph(MxVert), NV, NS, IPRINT
real(kind=wp), intent(in) :: XE(*), YE(*), ZE(*), RE(*), PTS(3,MxVert), CCC(3,MxVert)
real(kind=wp), intent(out) :: PP(3), AREA
integer(kind=iwp) :: I, JJ, N, N0, N1, N2, NSFE1
real(kind=wp) :: BETAN, COSPHIN, COSTN, DNORM, DNORM1, DNORM2, DNORM3, P1(3), P2(3), P3(3), PHIN, RR2, SCAL, SUM1, SUM2, U1(3), &
                 U2(3)
real(kind=wp), parameter :: Small = 1.0e-35_wp, TPI = Two*Pi

! Sfrutta il teorema di Gauss-Bonnet per calcolare l'area
! della tessera con vertici PTS(3,NV). Consideriamo sempre
! che il lato N della tessera e' quello compreso tra i vertici
! N e N+1 (oppure NV e 1). In CCC(3,NV) sono le posizioni dei
! centri degli archi che sottendono i vari lati della tessera.
! La formula di Gauss-Bonet per le sfere e':
!        Area=R^2[2pi+S(Phi(N)cosT(N))-S(Beta(N))]
! dove Phi(N) e' la lunghezza d'arco (in radianti) del lato N,
! T(N) e' l'angolo polare del lato N, Beta(N) l'angolo esterno
! relativo al vertice N.
!
! Questa SUBROUTINE viene chiamata anche da VOLCHG (in calcoli con
! dielettrico non omogeneo) per calcolare l'area delle tessere
! proiettate sulle relative "sfere del centro di carica": in questo caso
! i vertici della tessera sono in PTS, i centri degli archi (tutti
! coincidenti con il centro di carica) sono in CCC e il punto
! rappresentativo e' in PP. La cosa piu' importante e' che NS < 0.

N0 = 0 ! dummy initialize
N2 = 0 ! dummy initialize

! Calcola la prima sommatoria
SUM1 = Zero
do N=1,NV
  U1(:) = PTS(:,N)-CCC(:,N)
  if (N < NV) then
    U2(:) = PTS(:,N+1)-CCC(:,N)
  else
    U2(:) = PTS(:,1)-CCC(:,N)
  end if
  DNORM1 = U1(1)*U1(1)+U1(2)*U1(2)+U1(3)*U1(3)
  DNORM2 = U2(1)*U2(1)+U2(2)*U2(2)+U2(3)*U2(3)
  SCAL = U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3)
  COSPHIN = SCAL/(sqrt(DNORM1*DNORM2))
  if (COSPHIN > One) COSPHIN = One
  PHIN = acos(COSPHIN)

  ! Normalmente NS>0, ma se GAUBON e' stata chiamata da VOLCHG NS<0

  if (NS > 0) then
    ! NSFE1 e' la sfera con cui la sfera NS si interseca (eventualmente)
    NSFE1 = INTSPH(N)
    U1(1) = XE(NSFE1)-XE(NS)
    U1(2) = YE(NSFE1)-YE(NS)
    U1(3) = ZE(NSFE1)-ZE(NS)
  else
    U1(:) = ZERO
  end if
  DNORM1 = sqrt(U1(1)*U1(1)+U1(2)*U1(2)+U1(3)*U1(3))
  if (DNORM1 == ZERO) DNORM1 = One
  if (NS > 0) then
    U2(1) = PTS(1,N)-XE(NS)
    U2(2) = PTS(2,N)-YE(NS)
    U2(3) = PTS(3,N)-ZE(NS)
  else
    U2(:) = PTS(:,N)-CCC(:,1)
  end if
  DNORM2 = sqrt(U2(1)*U2(1)+U2(2)*U2(2)+U2(3)*U2(3))
  COSTN = (U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))/(DNORM1*DNORM2)
  SUM1 = SUM1+PHIN*COSTN
end do

! Calcola la seconda sommatoria: l'angolo esterno Beta(N) e'
! definito usando i versori (u(N-1),u(N)) tangenti alla sfera
! nel vertice N lungo le direzioni dei lati N-1 e N:
!            cos( Pi-Beta(N) )=u(N-1)*u(N)
!        u(N-1) = [V(N) x (V(N) x V(N-1))]/NORM
!        u(N) = [V(N) x (V(N) x V(N+1))]/NORM
! dove V(I) e' il vettore posizione del vertice I RISPETTO AL
! CENTRO DELL'ARCO CHE SI STA CONSIDERANDO.

SUM2 = Zero
! Loop sui vertici
do N=1,NV
  P1(:) = Zero
  N1 = N
  if (N > 1) N0 = N-1
  if (N == 1) N0 = NV
  if (N < NV) N2 = N+1
  if (N == NV) N2 = 1
  ! Trova i vettori posizione rispetto ai centri corrispondenti
  ! e i versori tangenti

  ! Lato N0-N1:
  P1(:) = PTS(:,N1)-CCC(:,N0)
  P2(:) = PTS(:,N0)-CCC(:,N0)

  call crprod(P1,P2,P3)
  P2(:) = P3(:)
  call crprod(P1,P2,P3)
  DNORM3 = sqrt(P3(1)*P3(1)+P3(2)*P3(2)+P3(3)*P3(3))
  if (DNorm3 < Small) DNorm3 = One
  U1(:) = P3(:)/DNORM3

  ! Lato N1-N2:
  P1(:) = PTS(:,N1)-CCC(:,N1)
  P2(:) = PTS(:,N2)-CCC(:,N1)

  call crprod(P1,P2,P3)
  P2(:) = P3(:)
  call crprod(P1,P2,P3)
  DNORM3 = sqrt(P3(1)*P3(1)+P3(2)*P3(2)+P3(3)*P3(3))
  if (DNorm3 < Small) DNorm3 = One
  U2(:) = P3(:)/DNORM3

  BETAN = acos(U1(1)*U2(1)+U1(2)*U2(2)+U1(3)*U2(3))
  SUM2 = SUM2+(PI-BETAN)
end do
! Calcola l'area della tessera
if (NS > 0) then
  AREA = RE(NS)*RE(NS)*(TPI+SUM1-SUM2)
else
  RR2 = (PTS(1,1)-CCC(1,1))**2+(PTS(2,1)-CCC(2,1))**2+(PTS(3,1)-CCC(3,1))**2
  AREA = RR2*(TPI+SUM1-SUM2)
end if
if (NS > 0) then
  ! Trova il punto rappresentativo (come media dei vertici)
  PP(:) = Zero
  do I=1,NV
    PP(1) = PP(1)+(PTS(1,I)-XE(NS))
    PP(2) = PP(2)+(PTS(2,I)-YE(NS))
    PP(3) = PP(3)+(PTS(3,I)-ZE(NS))
  end do
  DNORM = Zero
  do JJ=1,3
    DNORM = DNORM+PP(JJ)*PP(JJ)
  end do
  PP(1) = XE(NS)+PP(1)*RE(NS)/sqrt(DNORM)
  PP(2) = YE(NS)+PP(2)*RE(NS)/sqrt(DNORM)
  PP(3) = ZE(NS)+PP(3)*RE(NS)/sqrt(DNORM)
end if

! A causa delle approssimazioni numeriche, l'area di alcune piccole
! tessere puo' risultare negativa, e viene in questo caso trascurata
if (AREA < Zero) then
  AREA = Zero
  if (IPRINT >= 1) write(u6,1000) NS
end if

return

1000 format(/,'ATTENTION: THE SURFACE OF A TESSERA IN SPHERE ',I3,' IS NEGLECTED')

end subroutine GauBon
