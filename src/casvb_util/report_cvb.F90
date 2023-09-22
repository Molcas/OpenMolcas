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

subroutine report_cvb(orbs,norb)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: norb
real(kind=wp) :: orbs(norb,norb)
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstackr_cvb

write(u6,'(/,a)') ' Orbital coefficients :'
write(u6,'(a)') ' ----------------------'
call mxprint_cvb(orbs,norb,norb,0)
write(u6,'(/,a)') ' Overlap between orbitals :'
write(u6,'(a)') ' --------------------------'

i1 = mstackr_cvb(norb*norb)
call mxattb_cvb(orbs,orbs,norb,norb,norb,work(i1))
call mxprint_cvb(work(i1),norb,norb,0)
call mfreer_cvb(i1)

return

end subroutine report_cvb
