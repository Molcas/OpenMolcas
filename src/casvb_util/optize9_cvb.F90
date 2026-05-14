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

subroutine optize9_cvb(fx1,nparm,ioptc)

use casvb_global, only: formChk1, formChk2, formChk3
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: fx1
integer(kind=iwp), intent(in) :: nparm
integer(kind=iwp), intent(out) :: ioptc
integer(kind=iwp) :: ierr, iparm, it
real(kind=wp) :: cn, dum(1), e1, e2, fx
real(kind=wp), allocatable :: dx(:), grad(:), hessdx(:)
real(kind=wp), parameter :: tenth = 0.1_wp
real(kind=wp), external :: ddot_, rand_cvb

call mma_allocate(dx,nparm,label='dx')
call mma_allocate(grad,nparm,label='grad')
call mma_allocate(hessdx,nparm,label='hessdx')

call grad_cvb(grad)

dum(1) = rand_cvb(0.777_wp)
do iparm=1,nparm
  dx(iparm) = rand_cvb(Zero)-Half
end do
ierr = 0
call nize_cvb(dx,1,dum,nparm,0,ierr)
hessdx(:) = dx(:)
call hess_cvb(hessdx)

write(u6,'(a)') ' Simple check of gradient and Hessian using a random update vector :'
e1 = ddot_(nparm,dx,1,grad,1)
e2 = ddot_(nparm,dx,1,hessdx,1)
write(u6,'(a)') ' '
write(u6,formChk1) ' First-order change  :',e1
write(u6,formChk1) ' Second-order change :',e2
write(u6,'(a)') ' '

write(u6,formChk2) 'Norm     ','DFX(act) ','DFX(pred)','Ratio    ','F2(act)'
cn = One
do it=1,10
  call fxdx_cvb(fx,.false.,dx)
  write(u6,formChk3) cn,fx-fx1,cn*e1+cn*cn*Half*e2,(fx-fx1)/(cn*e1+cn*cn*Half*e2),(fx-fx1-cn*e1)/(cn*cn*Half)
  dx(:) = tenth*dx(:)
  cn = tenth*cn
end do

call mma_deallocate(dx)
call mma_deallocate(grad)
call mma_deallocate(hessdx)

ioptc = 0

return

end subroutine optize9_cvb
