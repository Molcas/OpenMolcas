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

subroutine asc2ab_cvb(detvec,nvec,nel,nalf)

use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: detvec(*)
integer(kind=iwp), intent(in) :: nvec, nel, nalf
integer(kind=iwp) :: nbet, ndet

call icomb_cvb(nel,nalf,ndet)
nbet = nel-nalf
call asc2ab2_cvb(detvec,nvec,nel,nalf,nbet,ndet)

return

end subroutine asc2ab_cvb
