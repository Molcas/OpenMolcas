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

function nPre(kS)

use Symmetry_Info, only: Mul
use input_mclr, only: nIsh, nOrb, nRS1, nRS2, nRS3, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nPre
integer(kind=iwp), intent(in) :: kS
integer(kind=iwp) :: iS, jS, nRest

nPre = 0
do iS=1,nSym
  jS = Mul(iS,kS)
  nRest = nOrb(jS)-nIsh(jS)
  nPre = nPre+nIsh(iS)*nRest**2
  nRest = nOrb(jS)-nRs1(jS)
  nPre = nPre+nRs1(iS)*nRest**2
  nRest = nOrb(jS)-nRs2(jS)
  nPre = nPre+nRs2(iS)*nRest**2
  nRest = nOrb(jS)-nRs3(jS)
  nPre = nPre+nRs3(iS)*nRest**2
end do

end function nPre
