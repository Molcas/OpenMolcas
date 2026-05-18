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

subroutine full2red(XLT,nXLT,Xab,nXab)

use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR, nSym
use constants, only: Zero
use definitions, only: iwp, wp

implicit none
integer(kind=iwp), intent(in) :: nXLT, nXab
real(kind=wp), intent(in) :: XLT(nXLT)
real(kind=wp), intent(out) :: Xab(nXab)
integer(kind=iwp) ISLT(8)
integer(kind=iwp), external :: cho_isao
integer(kind=iwp) iLoc, jSym, IS, ISYM, NB
integer(kind=iwp) jRab, kRab, iRab, iag, ibg, iSyma, ias, ibs, iab, kfrom

Xab(:) = Zero
! Select table column for use with caspt2:
iLoc = 3
! jSym=1 always: Used for density matrices.
jSym = 1
! Offsets to symmetry block in the LT matrix
IS = 0
do ISYM=1,NSYM
  ISLT(ISYM) = IS
  NB = NBAS(ISYM)
  IS = IS+(NB*(NB+1))/2
end do

do jRab=1,nnBstR(jSym,iLoc)
  kRab = iiBstr(jSym,iLoc)+jRab
  iRab = IndRed(kRab,iLoc)
  iag = iRS2F(1,iRab)
  ibg = iRS2F(2,iRab)
  iSyma = cho_isao(iag)
  ias = iag-ibas(iSyma)
  ibs = ibg-ibas(iSyma)
  if (ias >= ibs) then
    iab = (ias*(ias-1))/2+ibs
  else
    iab = (ibs*(ibs-1))/2+ias
  end if
  kfrom = isLT(iSyma)+iab
  Xab(jRab) = Xab(jRab)+XLT(kfrom)
end do

end subroutine full2red
