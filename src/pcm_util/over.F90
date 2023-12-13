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

subroutine Over(NSJ,ICOORD,GeoGrd,NAt,NTs,NEsf,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerRad,DerCentr)

use Constants, only: Zero, One, Two, Four, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSJ, ICOORD, NAt, NTs, NEsf, ISphe(NTs), NOrd(NEsf)
real(kind=wp), intent(out) :: GeoGrd
real(kind=wp), intent(in) :: Eps, Sphere(4,NEsf), Tessera(4,NTs), Q(2,NTs), DerRad(NEsf,NAt,3), DerCentr(NEsf,NAt,3,3)
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
