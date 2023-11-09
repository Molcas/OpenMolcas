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

subroutine clearcnt_cvb(icode)
! ICODE=1 : Orbitals changed
! ICODE=2 : CI coefficients changed
! ICODE=3 : Everything changed

use casvb_global, only: icnt_ci
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: icode
integer(kind=iwp) :: ichg, ipow1, ipow2

if (icode == 3) then
  icnt_ci(:) = 0
else
  ipow1 = 2
  ipow2 = 1
  do ichg=1,2
    if (mod(icode,ipow1) >= ipow2) icnt_ci(2:) = 0
    ipow1 = 2*ipow1
    ipow2 = 2*ipow2
  end do
end if

return

end subroutine clearcnt_cvb
