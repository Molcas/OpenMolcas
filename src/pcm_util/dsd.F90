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

subroutine dSd(IOpt,ITS,ICOORD,NESFJ,DERS,DERP,Tessera,Vert,Centr,nTs,Sphere,ISphe,IntSph,NewSph,NVert)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: MxVert = 20
integer(kind=iwp), intent(in) :: IOpt, ITS, ICOORD, NESFJ, nTs, ISphe(*), IntSph(MxVert,*), NewSph(2,*), NVert(*)
real(kind=wp), intent(out) :: DERS, DERP(3)
real(kind=wp), intent(in) :: Tessera(4,*), Vert(3,MxVert,*), Centr(3,MxVert,*), Sphere(4,*)
integer(kind=iwp) :: JJ, L, L0, L1, L2, L3, NSI, NV
real(kind=wp) :: DA, DERTS, DP(MxVert,3), PROD, PT(3), RSUM, SUMDER(3), SUMVERT(3)

! Trova le derivate dell'area e del punto rappresentativo della
! tessera ITS rispetto
!
!   IOpt = 0 : alla coordinata ICOORD(1=X, 2=Y, 3=Z) dell'atomo NSJ.
!              (ex dsdx)
!   IOpt = 1 : al raggio della sfera NSJ. (ex dsdr)
!
! NESFJ e' la sfera attorno all'atomo NSJ.
! ITS appartiene alla sfera NSI (diversa da NESFJ, e intersecata
! da NESFJ).

DERTS = ZERO
NSI = ISPHE(ITS)
NV = NVERT(ITS)
do L=1,NV
  DP(L,1) = ZERO
  DP(L,2) = ZERO
  DP(L,3) = ZERO
end do
do L=1,NV
  if (INTSPH(L,ITs) /= NESFJ) cycle
  ! L1 e' il lato di ITS che sta sulla sfera NESFJ
  L1 = L
  L0 = L-1
  if (L1 == 1) L0 = NV
  L2 = L1+1
  if (L1 == NV) L2 = 1
  L3 = L2+1
  if (L2 == NV) L3 = 1
  ! Trova la derivata del vertice L1: DP(L1)
  call DVer(IOpt,ICOORD,ITS,L0,L1,L2,DP(L1,1),DP(L1,2),DP(L1,3),Vert,Centr,nTs,Sphere,IntSph)
  ! Trova la derivata del vertice L2: DP(L2)
  call DVer(IOpt,ICOORD,ITS,L1,-L2,L3,DP(L2,1),DP(L2,2),DP(L2,3),Vert,Centr,nTs,Sphere,IntSph)

  ! 1) Derivata dell'area della tessera.

  ! Calcola i contributi dovuti ai lati L0-L1, L1-L2, L2-L3
  if (INTSPH(L0,ITs) /= NSI) then
    call DerPhi(IOpt,ICOORD,NESFJ,ITS,L0,L1,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
    DERTS = DERTS+DA
  end if
  call DerPhi(IOpt,ICOORD,NESFJ,ITS,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS+DA
  if (INTSPH(L2,ITs) /= ISPHE(ITS)) then
    call DerPhi(IOpt,ICOORD,NESFJ,ITS,L2,L3,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
    DERTS = DERTS+DA
  end if

  ! Calcola i contributi dovuti ai vertici L1 e L2
  call DerBet(IOpt,ICOORD,NESFJ,ITS,L0,L1,L2,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS-DA
  call DerBet(IOpt,ICOORD,NESFJ,ITS,L1,L2,L3,DP,DA,Vert,Centr,nTs,Sphere,IntSph,ISphe)
  DERTS = DERTS-DA
end do
DERS = DERTS

! 2) Derivata del punto rappresentativo.

! Trova il punto rappresentativo riferito al centro della sfera
PT(1) = Tessera(1,ITS)-Sphere(1,NSI)
PT(2) = Tessera(2,ITS)-Sphere(2,NSI)
PT(3) = Tessera(3,ITS)-Sphere(3,NSI)
! Trova il modulo della somma dei vertici della tessera
do JJ=1,3
  SUMVERT(JJ) = ZERO
end do
do L=1,NV
  SUMVERT(1) = SUMVERT(1)+(VERT(1,L,ITs)-Sphere(1,NSI))
  SUMVERT(2) = SUMVERT(2)+(VERT(2,L,ITs)-Sphere(2,NSI))
  SUMVERT(3) = SUMVERT(3)+(VERT(3,L,ITs)-Sphere(3,NSI))
end do
RSUM = ZERO
do JJ=1,3
  RSUM = RSUM+SUMVERT(JJ)*SUMVERT(JJ)
end do
RSUM = sqrt(RSUM)
! Trova la somma delle derivate dei vertici
do JJ=1,3
  SUMDER(JJ) = DP(L1,JJ)+DP(L2,JJ)
end do

PROD = ZERO
do JJ=1,3
  PROD = PROD+PT(JJ)*SUMDER(JJ)
end do
do JJ=1,3
  DERP(JJ) = (SUMDER(JJ)*Sphere(4,NSI)/RSUM)-(PT(JJ)*PROD/(Sphere(4,NSI)*RSUM))
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(NewSph)

end subroutine dSd
