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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mol2vbma_cvb(vecvb,vecmol,isyml,fac)

use Definitions, only: wp, iwp

implicit none
#include "main_cvb.fh"
real(kind=wp) :: vecvb(ndet), vecmol(*), fac
integer(kind=iwp) :: isyml
integer(kind=iwp) :: iwr, nsa, nsb
integer(kind=iwp), external :: mstackr_cvb

iwr = 2
call icomb_cvb(norb,nalf,nsa)
call icomb_cvb(norb,nbet,nsb)
call mol2vb2_cvb(vecvb,vecmol,isyml,fac,iwr,nsa,nsb)
ibasemx = max(ibasemx,mstackr_cvb(0))

return

end subroutine mol2vbma_cvb
