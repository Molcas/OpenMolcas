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

subroutine DerBet(IOpt,IC,NSFE,ITS,L0,L1,L2,DP,DA,Vert,Centr,Sphere,IntSph,ISphe)

use Constants, only: Zero, Two, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), parameter :: MxVert = 20
integer(kind=iwp), intent(in) :: IOpt, IC, NSFE, ITS, L0, L1, L2, IntSph(MxVert,*), ISphe(*)
real(kind=wp), intent(in) :: DP(MxVert,3), Vert(3,MxVert,*), Centr(3,MxVert,*), Sphere(4,*)
real(kind=wp), intent(out) :: DA
integer(kind=iwp) :: IDX, JJ, NS1, NS2, NS3
real(kind=wp) :: BETA, COSBETA, DBETA, DIST2, DNORM3, dnxxt0, dnxxt1, DT0(3), DT1(3), DV1(3), DV2(3), DV3(3), DV4(3), FACT, P1(3), &
                 P2(3), P3(3), PROD, SENBETA, T0(3), T1(3), T12(3), V1(3), V2(3), V3(3), V4(3)

! Trova la derivata dell'angolo esterno del vertice L1 della tessera
! ITS, formato dai lati L0-L1 e L1-L2,
!
! IOpt = 0 : rispetto alla coordinata IC (ex derbeta)
! IOpt = 1 : rispetto al raggio (ex derbeta1)
!
! della sfera NSFE, conoscendo la derivata della posizione del vertice P1.

NS1 = ISPHE(ITS)
NS2 = INTSPH(L0,ITs)
NS3 = INTSPH(L1,ITs)

! Trova i vettori posizione dei vertici L0, L1 rispetto al centro
! dell'arco L0, e dei vertici L1, L2 rispetto al centro dell'arco L1
do JJ=1,3
  V1(JJ) = VERT(JJ,L0,ITs)-CENTR(JJ,L0,ITs)
  V2(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L0,ITs)
  V3(JJ) = VERT(JJ,L1,ITs)-CENTR(JJ,L1,ITs)
  V4(JJ) = VERT(JJ,L2,ITs)-CENTR(JJ,L1,ITs)
end do
! Trova la derivata dei vettori posizione
do JJ=1,3
  DV1(JJ) = DP(L0,JJ)
  DV2(JJ) = DP(L1,JJ)
  DV3(JJ) = DP(L1,JJ)
  DV4(JJ) = DP(L2,JJ)
end do

! Corregge la derivata dei vettori posizione se il lato considerato
! appartiene alla sfera che si muove (distinguendo se ITS appartiene
! o no alla sfera che si muove)
IDX = 0
if ((NS1 == NSFE) .and. (NS2 /= NSFE)) IDX = 1
if ((NS1 /= NSFE) .and. (NS2 == NSFE)) IDX = 1
if (IDX == 1) then
  T12(1) = Sphere(1,NS2)-Sphere(1,NS1)
  T12(2) = Sphere(2,NS2)-Sphere(2,NS1)
  T12(3) = Sphere(3,NS2)-Sphere(3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    do JJ=1,3
      DV1(JJ) = DV1(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS2)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
      DV2(JJ) = DV2(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS2)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
    end do
    FACT = (Sphere(4,NS1)**2-Sphere(4,NS2)**2+DIST2)/(Two*DIST2)
    DV1(IC) = DV1(IC)-FACT
    DV2(IC) = DV2(IC)-FACT
  else if (IOpt == 1) then
    do JJ=1,3
      DV1(JJ) = DV1(JJ)+Sphere(4,NS2)*T12(JJ)/DIST2
      DV2(JJ) = DV2(JJ)+Sphere(4,NS2)*T12(JJ)/DIST2
    end do
  else
    write(u6,'(a)') 'Illegal IOpt in DerBet.'
    call Abend()
  end if
end if

IDX = 0
if ((NS1 == NSFE) .and. (NS3 /= NSFE)) IDX = 1
if ((NS1 /= NSFE) .and. (NS3 == NSFE)) IDX = 1
if (IDX == 1) then
  T12(1) = Sphere(1,NS3)-Sphere(1,NS1)
  T12(2) = Sphere(2,NS3)-Sphere(2,NS1)
  T12(3) = Sphere(3,NS3)-Sphere(3,NS1)
  DIST2 = T12(1)*T12(1)+T12(2)*T12(2)+T12(3)*T12(3)
  if (IOpt == 0) then
    do JJ=1,3
      DV3(JJ) = DV3(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS3)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
      DV4(JJ) = DV4(JJ)+(Sphere(4,NS1)**2-Sphere(4,NS3)**2)*T12(IC)*T12(JJ)/(DIST2*DIST2)
    end do
    FACT = (Sphere(4,NS1)**2-Sphere(4,NS3)**2+DIST2)/(Two*DIST2)
    DV3(IC) = DV3(IC)-FACT
    DV4(IC) = DV4(IC)-FACT
  else if (IOpt == 1) then
    do JJ=1,3
      DV3(JJ) = DV3(JJ)+Sphere(4,NS3)*T12(JJ)/DIST2
      DV4(JJ) = DV4(JJ)+Sphere(4,NS3)*T12(JJ)/DIST2
    end do
  else
    write(u6,'(a)') 'Illegal IOpt in DerBet.'
    call Abend()
  end if
end if

! Trova il vettore tangente al lato L0-L1: T0 = V2 x (V2 x V1)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  T0(JJ) = P3(JJ)
end do
dnxxt0 = DNORM3
! Trova il vettore tangente al lato L1-L2: T1 = V3 x (V3 x V4)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  T1(JJ) = P3(JJ)
end do
dnxxt1 = DNORM3

! Trova l'angolo Beta(L1)
PROD = T0(1)*T1(1)+T0(2)*T1(2)+T0(3)*T1(3)
BETA = PI-acos(PROD/(dnxxt0*dnxxt1))
COSBETA = cos(BETA)
SENBETA = sin(BETA)

! Trova la derivata della tangente T0 :
! DT0 = DV2 x (V2 x V1) + V2 x (DV2 x V1) + V2 x (V2 x DV1)

do JJ=1,3
  DT0(JJ) = Zero
end do

do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = DV2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = DV2(JJ)
  P2(JJ) = V1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = DV1(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V2(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT0(JJ) = DT0(JJ)+P3(JJ)
end do

! Trova la derivata della tangente T1 :
! DT1 = DV3 x (V3 x V4) + V3 x (DV3 x V4) + V3 x (V3 x DV4)

do JJ=1,3
  DT1(JJ) = Zero
end do

do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = DV3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = DV3(JJ)
  P2(JJ) = V4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = DV4(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  P1(JJ) = V3(JJ)
  P2(JJ) = P3(JJ)
end do
call VECP(P1,P2,P3,DNORM3)
do JJ=1,3
  DT1(JJ) = DT1(JJ)+P3(JJ)
end do

! Infine calcola la derivata dell'angolo Beta(L1)
do JJ=1,3
  P1(JJ) = T1(JJ)+COSBETA*dnxxt1*T0(JJ)/dnxxt0
  P2(JJ) = T0(JJ)+COSBETA*dnxxt0*T1(JJ)/dnxxt1
end do
do JJ=1,3
  DBETA = Zero
end do
do JJ=1,3
  DBETA = DBETA+(P1(JJ)*DT0(JJ)+P2(JJ)*DT1(JJ))
end do
DBETA = DBETA/(SENBETA*dnxxt0*dnxxt1)
DA = Sphere(4,NS1)*Sphere(4,NS1)*DBETA

return

end subroutine DerBet
