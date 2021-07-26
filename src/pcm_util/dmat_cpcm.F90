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

subroutine DMat_CPCM(iAt,iC,nTs,nS,nAt,fact,Tessera,DerMat,DerTes,DerPunt,DerCentr,iSphe)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAt, iC, nTs, nS, nAt, iSphe(nTs)
real(kind=wp), intent(in) :: fact, Tessera(4,nTs), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), DerCentr(nS,nAt,3,3)
real(kind=wp), intent(out) :: DerMat(nTs,nTs)
integer(kind=iwp) :: ITs, JTs, L, LJ
real(kind=wp) :: DIJ, DXIJ, DYIJ, DZIJ, PROD, XIJ, YIJ, ZIJ

! Compute the derivative of the CPCM matrix wrt atom iat, coord. ic

! Loop on tesserae
do ITs=1,NTs
  L = ISPHE(ITs)
  do JTs=1,NTs
    LJ = ISPHE(JTs)
    ! Diagonal elements
    if (ITs == JTs) then
      DerMat(ITs,ITs) = fact*DERTES(ITs,iAt,IC)/(Tessera(4,ITs)*sqrt(Tessera(4,ITs)))
    else
      ! Off diagonal elements
      XIJ = Tessera(1,ITs)-Tessera(1,JTs)
      YIJ = Tessera(2,ITs)-Tessera(2,JTs)
      ZIJ = Tessera(3,ITs)-Tessera(3,JTs)
      DIJ = sqrt(XIJ*XIJ+YIJ*YIJ+ZIJ*ZIJ)
      DXIJ = DERPUNT(ITs,iAt,IC,1)+DERCENTR(L,iAt,IC,1)-DERPUNT(JTs,iAt,IC,1)-DERCENTR(LJ,iAt,IC,1)
      DYIJ = DERPUNT(ITs,iAt,IC,2)+DERCENTR(L,iAt,IC,2)-DERPUNT(JTs,iAt,IC,2)-DERCENTR(LJ,iAt,IC,2)
      DZIJ = DERPUNT(ITs,iAt,IC,3)+DERCENTR(L,iAt,IC,3)-DERPUNT(JTs,iAt,IC,3)-DERCENTR(LJ,iAt,IC,3)
      PROD = (XIJ*DXIJ+YIJ*DYIJ+ZIJ*DZIJ)/DIJ**3
      DerMat(ITs,JTs) = -PROD
    end if
  end do
end do

return

end subroutine DMat_CPCM
