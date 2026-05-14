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

subroutine mol2vb_cvb(vecvb,vecmol,isyml)

use casvb_global, only: nalf, nbet, ndet, norb
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: vecvb(ndet)
real(kind=wp), intent(_IN_) :: vecmol(*)
integer(kind=iwp), intent(in) :: isyml
integer(kind=iwp) :: iwr, nsa, nsb

iwr = 1
call icomb_cvb(norb,nalf,nsa)
call icomb_cvb(norb,nbet,nsb)
call mol2vb2_cvb(vecvb,vecmol,isyml,Zero,iwr,nsa,nsb)

return

end subroutine mol2vb_cvb
