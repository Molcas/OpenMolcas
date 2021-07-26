!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Giovanni Ghigo                                         *
!***********************************************************************

subroutine LenInt(iSymI,iSymJ,iSymA,iSymB,nProdIJ,nProdAB,nProdE1,nProdE2)
!***********************************************************************
! Author  :  Giovanni Ghigo                                            *
!            Lund University, Sweden                                   *
!----------------------------------------------------------------------*
! Return the Length of Coulomb (nProdAB), Exchanges (nProdE1 and       *
! nProdE2) matrices for each i,j and the length (nProdIJ) of the i,j   *
! matrix for each Symmetry Block (iSymI,iSymJ,iSymA,iSymB)             *
!***********************************************************************

use Cho_Tra, only: DoTCVA, nOrb, nOsh, nSsh
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iSymI, iSymJ, iSymA, iSymB
integer(kind=iwp), intent(out) :: nProdIJ, nProdAB, nProdE1, nProdE2
integer(kind=iwp) :: nExtA, nExtB, nOccI, nOccJ, nOrbA, nOrbB

nProdIJ = 0
nProdAB = 0
nProdE1 = 0
nProdE2 = 0
nOccI = nOsh(iSymI)
nOccJ = nOsh(iSymJ)
nOrbA = nOrb(iSymA)
nOrbB = nOrb(iSymB)
nExtA = nSsh(iSymA)
nExtB = nSsh(iSymB)
if (iSymI == iSymJ) then
  nProdIJ = nOccI*(nOccJ+1)/2
else
  nProdIJ = nOccI*nOccJ
end if
if (iSymA == iSymB) then
  nProdAB = nOrbA*(nOrbB+1)/2
else
  if (iSymA > iSymB) then
    nProdAB = nOrbA*nOrbB
  else
    nProdAB = 0
  end if
end if
if (iSymA >= iSymB) then
  if (DoTCVA) then
    nProdE1 = nOrbA*nOrbB
  else
    nProdE1 = nExtA*nExtB
  end if
  nProdE2 = 0
else
  nProdE1 = 0
  if (DoTCVA) then
    nProdE2 = nOrbA*nOrbB
  else
    nProdE2 = nExtA*nExtB
  end if
end if

return

end subroutine LenInt
