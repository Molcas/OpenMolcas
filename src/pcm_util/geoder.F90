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

subroutine GeoDer(nAt,Cond,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerDM,Grd,DerTes,DerPunt,DerRad,DerCentr)

use Constants, only: Zero, Half, Angstrom
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, nS, ISphe(*), NOrd(*)
logical(kind=iwp), intent(in) :: Cond
real(kind=wp), intent(in) :: Eps, Sphere(4,*), Tessera(4,*), Q(2,*), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), DerRad(nS,nAt,3), &
                             DerCentr(nS,nAt,3,3)
real(kind=wp), intent(_OUT_) :: DerDM(nTs,*), Grd(3,*)
integer(kind=iwp) :: IAtom, iTs, IXYZ, jTs
real(kind=wp) :: GeoGrd, Qi, Qj

! Compute the PCM geometric contribution to gradients

call FZero(Grd,3*nAt)
call FZero(DerDM,nTs*nTs)
do IAtom=1,nAt
  do IXYZ=1,3
    ! Dielectric model
    if (.not. Cond) then
      call Over(IAtom,IXYZ,GeoGrd,nAt,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerRad,DerCentr)
    ! Conductor model
    elseif (Cond) then
      GeoGrd = Zero
      call DerD(Angstrom,IAtom,IXYZ,Tessera,ISphe,DerDM,DerTes,DerPunt,DerCentr,nTs,nAt,nS)
      do iTs=1,nTs
        Qi = Q(1,iTs)+Q(2,iTs)
        do jTs=1,nTs
          Qj = Q(1,jTs)+Q(2,jTs)
          GeoGrd = GeoGrd+Qi*DerDM(iTs,jTs)*Qj
        end do
      end do
    end if
    Grd(IXYZ,IAtom) = GeoGrd*Half
  end do
end do

return

end subroutine GeoDer
!====
subroutine DERD(ToAng,NSJ,ICOORD,Tessera,ISphe,DerDM,DerTes,DerPunt,DerCentr,NTs,NAt,NEsf)

use Constants, only: One, Pi
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSJ, ICOORD, ISphe(*), NTs, NAt, NEsf
real(kind=wp), intent(in) :: ToAng, Tessera(4,*), DerTes(NTs,NAt,*), DerPunt(NTs,NAt,3,3), DerCentr(NEsf,NAt,3,3)
real(kind=wp), intent(_OUT_) :: DerDM(NTs,*)
integer(kind=iwp) :: ITS, JTS, L, LJ
real(kind=wp) :: ANTOAU, DIJ, DXIJ, DYIJ, DZIJ, FAC, PROD, XIJ, YIJ, ZIJ

ANTOAU = One/ToAng

! Calcola la matrice DERDM, derivata di DMAT (secondo COSMO)
! rispetto alla coordinata ICOORD della sfera NSJ

! Loop sugli elementi di DMAT
do ITS=1,NTS
  L = ISPHE(ITS)
  do JTS=1,NTS
    LJ = ISPHE(JTS)
    ! Elementi diagonali di DMAT(x)
    if (ITS == JTS) then
      FAC = -1.0694_wp*sqrt(PI)
      DERDM(ITS,ITS) = FAC*DERTES(ITS,NSJ,ICOORD)*ANTOAU/(Tessera(4,ITS)*sqrt(Tessera(4,ITS)))
    else
      ! Elementi fuori diagonale di DMAT(x)
      XIJ = Tessera(1,ITS)-Tessera(1,JTS)
      YIJ = Tessera(2,ITS)-Tessera(2,JTS)
      ZIJ = Tessera(3,ITS)-Tessera(3,JTS)
      DIJ = sqrt(XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ)
      DXIJ = DERPUNT(ITS,NSJ,ICOORD,1)+DERCENTR(L,NSJ,ICOORD,1)-DERPUNT(JTS,NSJ,ICOORD,1)-DERCENTR(LJ,NSJ,ICOORD,1)
      DYIJ = DERPUNT(ITS,NSJ,ICOORD,2)+DERCENTR(L,NSJ,ICOORD,2)-DERPUNT(JTS,NSJ,ICOORD,2)-DERCENTR(LJ,NSJ,ICOORD,2)
      DZIJ = DERPUNT(ITS,NSJ,ICOORD,3)+DERCENTR(L,NSJ,ICOORD,3)-DERPUNT(JTS,NSJ,ICOORD,3)-DERCENTR(LJ,NSJ,ICOORD,3)
      PROD = (XIJ*DXIJ+YIJ*DYIJ+ZIJ*DZIJ)/DIJ**3
      DERDM(ITS,JTS) = -PROD
    end if
  end do
end do

return

end subroutine DERD
!====
subroutine Over(NSJ,ICOORD,GeoGrd,NAt,NTs,NEsf,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerRad,DerCentr)

use Constants, only: Zero, One, Two, Four, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSJ, ICOORD, NAt, NTs, NEsf, ISphe(*), NOrd(*)
real(kind=wp), intent(out) :: GeoGrd
real(kind=wp), intent(in) :: Eps, Sphere(4,*), Tessera(4,*), Q(2,*), DerRad(NEsf,NAt,3), DerCentr(NEsf,NAt,3,3)
integer(kind=iwp) :: I, ITS, L, NESFJ
real(kind=wp) :: DCENTN, DN, Fact, SESE, SESN, SNSN, XNI, YNI, ZNI

! Calcola la sovrapposizione delle densita' superficiali sulla
! porzione di superficie che appartiene alla sfera
! contenente l'atomo NSJ rispetto al quale stiamo derivando.
!
! NESFJ e' la sfera che sta attorno all'atomo NSJ: se NSJ non ha
! nessuna sfera, NESFJ = 0

NESFJ = 0
do I=1,NESF
  if (NSJ == NORD(I)) NESFJ = I
end do

SESE = ZERO
SNSN = ZERO
SESN = ZERO
do ITS=1,NTS
  L = ISPHE(ITS)
  XNI = -(Sphere(1,L)-Tessera(1,ITS))/Sphere(4,L)
  YNI = -(Sphere(2,L)-Tessera(2,ITS))/Sphere(4,L)
  ZNI = -(Sphere(3,L)-Tessera(3,ITS))/Sphere(4,L)

  if (L == NESFJ) then
    DN = ZERO ! dummy initialize
    if (ICOORD == 1) DN = XNI
    if (ICOORD == 2) DN = YNI
    if (ICOORD == 3) DN = ZNI
  else
    DCENTN = XNI*DerCentr(L,NSJ,ICOORD,1)+YNI*DerCentr(L,NSJ,ICOORD,2)+ZNI*DerCentr(L,NSJ,ICOORD,3)
    DN = DERRAD(L,NSJ,ICOORD)+DCENTN
  end if

  SNSN = SNSN+DN*Q(1,ITS)**2/Tessera(4,ITS)
  SESE = SESE+DN*Q(2,ITS)**2/Tessera(4,ITS)
  SESN = SESN+DN*Q(1,ITS)*Q(2,ITS)/Tessera(4,ITS)
end do
Fact = Four*PI*Eps/(Eps-One)
GeoGrd = Fact*(SESE+SNSN+Two*SESN)

return

end subroutine Over
