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

subroutine Der_Norm(iAt1,iCoord1,iAt2,iCoord2,nTs,nAt,nS,Tessera,Der1,DerRad,DerTes,DerPunt,Sphere,iSphe,nOrd)

use Constants, only: Zero, One, Angstrom
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iAt1, iCoord1, iAt2, iCoord2, nTs, nAt, nS, iSphe(nTs), nOrd(nS)
real(kind=wp), intent(in) :: Tessera(4,nTs), DerRad(nS,nAt,3), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), Sphere(4,nS)
real(kind=wp), intent(out) :: Der1(nTs)
integer(kind=iwp) :: iAt1_S, iS, iTs, L
real(kind=wp) :: AnToAU, Der_1, dN, dN_Its, SqArea

AnToAU = One/Angstrom
! Find out if atom iAt1 has a sphere around
iAt1_S = 0
do iS=1,nS
  if (iAt1 == nOrd(iS)) iAt1_S = iS
end do

! Compute the derivative of the normal vector over the tessera area

! Loop on tiles
dN = Zero ! Dummy initialization.
Der1(:) = Zero
do iTs=1,nTs
  L = iSphe(iTs)
  if (L == iAt1_S) then
    if (iCoord1 == 1) dN = (Sphere(1,L)-Tessera(1,iTs))/Sphere(4,L)
    if (iCoord1 == 2) dN = (Sphere(2,L)-Tessera(2,iTs))/Sphere(4,L)
    if (iCoord1 == 3) dN = (Sphere(3,L)-Tessera(3,iTs))/Sphere(4,L)
    dN_Its = DerPunt(iTs,iAt2,iCoord2,iCoord1)+dN*DerRad(L,iAt2,iCoord2)
    dN_Its = -dN_Its/Sphere(4,L)
  else
    dN = Zero
    dN_Its = Zero
  end if
  SqArea = Tessera(4,iTs)*Tessera(4,iTs)
  Der_1 = -dN*DerTes(iTs,iAt2,iCoord2)*AnToAU/SqArea
  Der1(iTs) = -dN_Its/Tessera(4,iTs)-Der_1
end do

return

end subroutine Der_Norm
