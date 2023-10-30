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

subroutine makecivbs_cvb(civbs,orbs,cvbdet)
! Construct CIVBS ( = T(s) * CIVB ):

use casvb_global, only: icnt_ci, ndet, ndetvb, norb
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: civbs(0:ndet), cvbdet(ndetvb)
real(kind=wp), intent(in) :: orbs(norb,norb)
integer(kind=iwp) :: icivbs

icivbs = nint(civbs(0))
if (icnt_ci(icivbs) == 4) return

call vb2cic_cvb(cvbdet,civbs)
call applyts_cvb(civbs,orbs)
icnt_ci(icivbs) = 4

return

end subroutine makecivbs_cvb
