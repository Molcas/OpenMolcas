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

subroutine FnoMP2_putInf(mSym,lnOrb,lnOcc,lnFro,lnDel,lnVir,ip_X,ip_Y)
! Purpose: put info in MP2 common blocks.

#include "implicit.fh"
integer lnOrb(8), lnOcc(8), lnFro(8), lnDel(8), lnVir(8)
integer ip_X, ip_Y
#include "corbinf.fh"
#include "chomp2_cfg.fh"

nSym = mSym

do iSym=1,nSym
  nOrb(iSym) = lnOrb(iSym)
  nOcc(iSym) = lnOcc(iSym)
  nFro(iSym) = lnFro(iSym)
  nDel(iSym) = lnDel(iSym)
  nExt(iSym) = lnVir(iSym)
end do

DoFNO = .true.
ip_Dab = ip_X
ip_Dii = ip_Y
l_Dab = nExt(1)
l_Dii = nOcc(1)
do iSym=2,nSym
  l_Dab = l_Dab+nExt(iSym)**2
  l_Dii = l_Dii+nOcc(iSym)
end do

return

end subroutine FnoMP2_putInf
