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
! Copyright (C) 2008, Francesco Aquilante                              *
!***********************************************************************

subroutine FnoMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
! Purpose: put info in MP2 common blocks.

use ChoMP2, only: DoFNO, l_Dii
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mSym, lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
integer(kind=iwp) :: iSym
#include "corbinf.fh"

nSym = mSym

do iSym=1,nSym
  nOrb(iSym) = lnOrb(iSym)
  nOcc(iSym) = lnOcc(iSym)
  nFro(iSym) = lnFro(iSym)
  nDel(iSym) = lnDel(iSym)
  nExt(iSym) = lnVir(iSym)
end do

DoFNO = .true.
l_Dii = nOcc(1)
do iSym=2,nSym
  l_Dii = l_Dii+nOcc(iSym)
end do

return

end subroutine FnoMP2_putInf
