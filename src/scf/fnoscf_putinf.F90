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
use cOrbInf, only: nDel, nExt, nFro, nOcc, nOrb, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: mSym, lnOrb(mSym), lnOcc(mSym), lnFro(mSym), lnDel(mSym), lnVir(mSym)

nSym = mSym

nOrb(1:nSym) = lnOrb(:)
nOcc(1:nSym) = lnOcc(:)
nFro(1:nSym) = lnFro(:)
nDel(1:nSym) = lnDel(:)
nExt(1:nSym) = lnVir(:)

DoFNO = .true.
l_Dii = sum(nOcc(:))

end subroutine FnoSCF_putInf
