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

subroutine mksymcvb_cvb()

use casvb_global, only: cvb, iconstruc, ipr, nvb, vbdet
use Definitions, only: wp, u6

implicit none
real(kind=wp) :: psnrm
real(kind=wp), parameter :: thresh = 1.0e-15_wp
real(kind=wp), external :: ddot_

! Constraints on struc coeffs - either symmetry or deleted
if (iconstruc > 0) then
  if (ipr(1) >= 0) write(u6,'(/,a)') ' Imposing constraints on the structure coefficients.'
  call symtrizcvb_cvb(cvb)
  psnrm = ddot_(nvb,cvb,1,cvb,1)
  if (psnrm < thresh) then
    write(u6,*) ' Fatal error - structure coefficients null after symmetrization!'
    call abend_cvb()
  end if
  if (ipr(1) >= 0) then
    write(u6,'(/,a)') ' Constrained structure coefficients :'
    write(u6,'(a)') ' ------------------------------------'
    call vecprint_cvb(cvb,nvb)
  end if
end if
call str2vbc_cvb(cvb,vbdet)

return

end subroutine mksymcvb_cvb
