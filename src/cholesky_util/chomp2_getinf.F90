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
! Copyright (C) 2004, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine ChoMP2_GetInf(lnOrb,lnOcc,lnFro,lnDel,lnVir)
!
! Thomas Bondo Pedersen, Dec. 2004.
!
! Purpose: get info from conventional MP2 common blocks.

implicit real*8(a-h,o-z)
integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
#include "corbinf.fh"

do iSym=1,nSym
  lnOrb(iSym) = nOrb(iSym)
  lnOcc(iSym) = nOcc(iSym)
  lnFro(iSym) = nFro(iSym)
  lnDel(iSym) = nDel(iSym)
  lnVir(iSym) = nExt(iSym)
end do

end subroutine ChoMP2_GetInf
