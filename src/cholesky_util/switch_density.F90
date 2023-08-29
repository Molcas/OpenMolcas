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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine switch_density(iLoc,XLT,Xab,kSym)

use Index_Functions, only: iTri
use Cholesky, only: iBas, iiBstR, IndRed, iRS2F, nnBstR
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: iLoc, kSym
real(kind=wp), intent(in) :: XLT(*)
real(kind=wp), intent(_OUT_) :: Xab(*)
integer(kind=iwp) :: iab, iag, ias, ibg, ibs, iRab, iSyma, jRab, jSym, kRab
real(kind=wp) :: xf
integer(kind=iwp), external :: cho_isao

jSym = 1 ! only total symmetric density

do jRab=1,nnBstR(jSym,iLoc)

  kRab = iiBstr(jSym,iLoc)+jRab
  iRab = IndRed(kRab,iLoc)

  iag = iRS2F(1,iRab)  !global address
  ibg = iRS2F(2,iRab)

  iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

  xf = real(1-min(1,abs(iSyma-kSym)),kind=wp)

  ias = iag-ibas(iSyma)  !address within that symm block
  ibs = ibg-ibas(iSyma)
  iab = iTri(ias,ibs)

  Xab(jRab) = xf*XLT(iab)

end do  ! jRab loop

return

end subroutine switch_density
