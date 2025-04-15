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
integer(kind=iwp) :: kS
integer(kind=iwp) :: iOut, iS, jS, nRest

iOut = 0
do is=1,nSym
  jS = Mul(iS,kS)
  nRest = nOrb(js)-nIsh(js)
  iOut = iOut+nIsh(is)*nRest*(nRest+1)
  nRest = nOrb(js)-nRs1(js)
  iOut = iOut+nRs1(is)*nRest*(nRest+1)
  nRest = nOrb(js)-nRs2(js)
  iOut = iOut+nRs2(is)*nRest*(nRest+1)
  nRest = nOrb(js)-nRs3(js)
  iOut = iOut+nRs3(is)*nRest*(nRest+1)
end do
nPre = iOut

end function nPre
