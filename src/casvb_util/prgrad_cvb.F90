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

subroutine prgrad_cvb(grad,n)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: n
real(kind=wp) :: grad(n)
#include "main_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i1
integer(kind=iwp), external :: mstackr_cvb

if (ip(3) < 2) return
i1 = mstackr_cvb(norb*norb)
call mxunfold_cvb(grad,work(i1),norb)
write(u6,'(/,a)') ' Orbital gradient :'
call mxprint_cvb(work(i1),norb,norb,0)
if (n-nprorb > 0) then
  write(u6,'(a)') ' Structure coefficient gradient :'
  call mxprint_cvb(grad(nprorb+1),1,n-nprorb,0)
end if
call mfreer_cvb(i1)

return

end subroutine prgrad_cvb
