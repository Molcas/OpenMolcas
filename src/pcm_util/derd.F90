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

subroutine DERD(NSJ,ICOORD,Tessera,ISphe,DerDM,DerTes,DerPunt,DerCentr,NTs,NAt,NEsf)

use PCM_Arrays, only: DiagScale
use Constants, only: One, Pi, Angstrom
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NSJ, ICOORD, NTs, ISphe(NTs), NAt, NEsf
real(kind=wp), intent(in) :: Tessera(4,NTs), DerTes(NTs,NAt,3), DerPunt(NTs,NAt,3,3), DerCentr(NEsf,NAt,3,3)
real(kind=wp), intent(out) :: DerDM(NTs,NTs)
integer(kind=iwp) :: ITS, JTS, L, LJ
real(kind=wp) :: ANTOAU, DIJ, DXIJ, DYIJ, DZIJ, FAC, PROD, XIJ, YIJ, ZIJ

ANTOAU = One/Angstrom

! Calcola la matrice DERDM, derivata di DMAT (secondo COSMO)
! rispetto alla coordinata ICOORD della sfera NSJ

! Loop sugli elementi di DMAT
do ITS=1,NTS
  L = ISPHE(ITS)
  do JTS=1,NTS
    LJ = ISPHE(JTS)
    ! Elementi diagonali di DMAT(x)
    if (ITS == JTS) then
      FAC = -DiagScale*sqrt(PI)
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
