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
! Copyright (C) 2010, Francesco Aquilante                              *
!***********************************************************************

subroutine FnoSCF_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir)
! Purpose: put info in MP2 common blocks.

use ChoMP2, only: DoFNO, l_Dii
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: mSym, lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
#include "corbinf.fh"

nSym = mSym

nOrb(1:nSym) = lnOrb(1:nSym)
nOcc(1:nSym) = lnOcc(1:nSym)
nFro(1:nSym) = lnFro(1:nSym)
nDel(1:nSym) = lnDel(1:nSym)
nExt(1:nSym) = lnVir(1:nSym)

DoFNO = .true.
l_Dii = sum(nOcc(1:nSym))

return

end subroutine FnoSCF_putInf
