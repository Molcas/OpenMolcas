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

subroutine gr_svb1_cvb(civecp,civbs,civb,dvbdet,grad,grad1,grad2,gradx,vec1)

use casvb_global, only: aa1, aa2, ndet, ndetvb, nfr, norb, npr, oaa2, oaa3, ovraa, ovrab
use Constants, only: Zero, One, Two, Three, Four
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: civecp(0:ndet), civbs(0:ndet)
real(kind=wp), intent(in) :: civb(0:ndet)
real(kind=wp), intent(out) :: dvbdet(ndetvb), grad(nfr), grad1(npr), grad2(npr), gradx(norb,norb), vec1(npr)

aa1 = One/sqrt(ovraa)
aa2 = -aa1/(Two*ovraa)
oaa2 = Two*ovrab*aa2
oaa3 = Three*ovrab*aa1/(Four*ovraa*ovraa)

gradx(:,:) = Zero
call onedens_cvb(civb,civbs,gradx,.true.,1)

call mkgrd_cvb(civb,civbs,grad1,dvbdet,npr,.true.)
call mkgrd_cvb(civb,civecp,grad2,dvbdet,npr,.true.)

vec1(:) = aa1*grad2(:)+oaa2*grad1(:)
grad1(:) = Two*grad1(:)
call prgrad_cvb(vec1,npr)

call make_cvb('ORBFREE')
call make_cvb('CIFREE')
call all2free_cvb(vec1,grad,1)

return

end subroutine gr_svb1_cvb
