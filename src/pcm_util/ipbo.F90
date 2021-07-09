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

function IPBO(ToAng,IA,JA,RIJ,BondOr)
! Return the order of the bond between atoms of atomic numbers IA
! and JA, separated by RIJ, or 0 if there is no bond.

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IPBO
real(kind=wp), intent(in) :: ToAng, RIJ
integer(kind=iwp), intent(in) :: IA, JA
real(kind=wp), intent(out) :: BondOr
integer(kind=iwp) :: IBondO
real(kind=wp) :: R0IJ, R1IJ
real(kind=wp), external :: RCov97

! Generate connectivity based on bond distances alone.  The criteria
! is whether the distances is no more than 30% longer than the sum of
! covalent radii of involved atoms. For the moment all bond types are
! determined using Pauling bond orders.
! Note that RIJ is multiplied by ToAng, i.e. if it is already in Ang.
! ToAng must be set = 1

IPBO = 0
R1IJ = RIJ*ToAng
R0IJ = RCov97(IA,JA)
!if (R1IJ > R0IJ*1.3_wp) return
BondOr = exp((R0IJ-R1IJ)/0.3_wp)
if (BondOr < 0.2_wp) return
IBondO = int(BondOr+Half)
IBondO = max(IBondO,1)
IBondO = min(IBondO,3)
IPBO = IBondO

return

end function IPBO
