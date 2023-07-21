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

use ChoArr, only: iRS2F
use ChoSwp, only: IndRed
use stdalloc

implicit real*8(a-h,o-z)
integer, external :: cho_isao
integer iLoc, kSym
real*8 XLT(*), Xab(*)
#include "cholesky.fh"
#include "choorb.fh"
! Statement function
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

jSym = 1 ! only total symmetric density

do jRab=1,nnBstR(jSym,iLoc)

  kRab = iiBstr(jSym,iLoc)+jRab
  iRab = IndRed(kRab,iLoc)

  iag = iRS2F(1,iRab)  !global address
  ibg = iRS2F(2,iRab)

  iSyma = cho_isao(iag)  !symmetry block; Sym(b)=Sym(a)

  xf = dble(1-min(1,abs(iSyma-kSym)))

  ias = iag-ibas(iSyma)  !address within that symm block
  ibs = ibg-ibas(iSyma)
  iab = iTri(ias,ibs)

  Xab(jRab) = xf*XLT(iab)

end do  ! jRab loop

return

end subroutine switch_density
