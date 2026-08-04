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

subroutine red2full(XLT,nXLT,Xab,nXab)

use Index_Functions, only: iTri, nTri_Elem
use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nBas, nnBstR, nSym
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nXLT, nXab
real(kind=wp), intent(inout) :: XLT(nXLT)
real(kind=wp), intent(in) :: Xab(nXab)
integer(kind=iwp) :: iab, iag, ias, ibg, ibs, iLoc, iRab, IS, ISLT(8), ISYM, iSyma, jRab, jSym, kRab, kto, NB
integer(kind=iwp), external :: cho_isao

! Select table column for use with caspt2:
iLoc = 3
! jSym=1 always: Used for density matrices.
jSym = 1
! Offsets to symmetry block in the LT matrix
IS = 0
do ISYM=1,NSYM
  ISLT(ISYM) = IS
  NB = NBAS(ISYM)
  IS = IS+nTri_Elem(NB)
end do

do jRab=1,nnBstR(jSym,iLoc)
  kRab = iiBstr(jSym,iLoc)+jRab
  iRab = IndRed(kRab,iLoc)
  iag = iRS2F(1,iRab)
  ibg = iRS2F(2,iRab)
  iSyma = cho_isao(iag)
  ias = iag-ibas(iSyma)
  ibs = ibg-ibas(iSyma)
  iab = iTri(ias,ibs)
  kto = isLT(iSyma)+iab
  XLT(kto) = XLT(kto)+Xab(jRab)
end do

end subroutine red2full
