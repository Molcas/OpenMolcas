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

subroutine gr_evb1_cvb(civbh,civbs,civb,dvbdet,grad,grad1,grad2,gradx,vec1)

use casvb_global, only: f1, f2, f3, f4, ndet, ndetvb, nfr, norb, npr, ovraa, ww
use Constants, only: Zero, One, Two
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civbh(0:ndet), civbs(0:ndet)
real(kind=wp), intent(in) :: civb(0:ndet)
real(kind=wp), intent(out) :: dvbdet(ndetvb), grad(nfr), grad1(npr), grad2(npr), gradx(norb,norb), vec1(npr)

f1 = One/ovraa
f2 = -Two*f1*f1*ww
f1 = Two*f1
f3 = -f1*f1
f4 = -f1*f3*ww

gradx(:,:) = Zero
call onedens_cvb(civb,civbs,gradx,.true.,1)

call mkgrd_cvb(civb,civbs,grad1,dvbdet,npr,.true.)
call mkgrd_cvb(civb,civbh,grad2,dvbdet,npr,.true.)

! Use VEC1 as work:
vec1(:) = f1*grad2(:)+f2*grad1(:)
call prgrad_cvb(vec1,npr)

call make_cvb('ORBFREE')
call make_cvb('CIFREE')
call all2free_cvb(vec1,grad,1)

return

end subroutine gr_evb1_cvb
