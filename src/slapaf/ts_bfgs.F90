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
! Copyright (C) 2005, Christian Ander                                  *
!***********************************************************************

subroutine ts_bfgs(dq,y,H,nH)
!     Hessian update method; TS-BFGS ; Bofill - "Remarks on the
!     Updated Hessian Matrix Methods" 2003.
!
!     Vectors in the article are column vectors.
!
!     Implemented by Christian Ander, 2005, christian@eximius.se
!                                                                      *
!***********************************************************************
!                                                                      *
!     dq    :  Perturbation Geometry      (nH)
!     y     :  Difference in Gradient     (nH)
!     WorkM :  Temporary working matrix   (nH,nH)
!     WorkV :  Temporary working vector   (nH)
!     WorkR :  Temporary working variable (real*8)
!     H     :  Hessian                    (nH,nH)
!     nH    :  Hessian size               (integer)
!     f,a,b :  Multi-used variables       (real*8)
!     v,u   :  Multi-used vectors         (nH)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nH
real(kind=wp), intent(in) :: dq(nH), y(nH)
real(kind=wp), intent(inout) :: H(nH,nH)
integer(kind=iwp) :: i
real(kind=wp) :: a, b, f, WorkR
real(kind=wp), allocatable :: u(:), v(:), WorkM(:,:), WorkV(:)
real(kind=wp), external :: ddot_

call mma_allocate(WorkM,nH,nH,label='WorkM')
call mma_allocate(WorkV,nH,label='WorkV')
call mma_allocate(v,nH,label='v')
call mma_allocate(u,nH,label='u')
!                                                                      *
!***********************************************************************
!                                                                      *
! Equation 23 calculates the new Hessian, some simplifications :
!
! H_k+1 = H_k + 1/f (v*u^T + u*v^T - (y^T*dq - dq^T*B*dq)/f * u*u^T)
!
! Where f = (y^T*dq)^2 + (dq^T|B|dq)^2 = a^2 + b^2
! ,     u = y^T*dq*y + (dq^T|B|dq)|B|dq = a*y + b|B|dq
! and   v = y - B*dq (quasi-Newton condition)

#ifdef _DEBUGPRINT_
! Make a comment in logfile
write(u6,*) 'hello from ts_bfgs'
#endif

! Calculation of u = y^T*dq*y + (dq^T|B|dq)|B|dq

! a = y^T * dq

a = DDot_(nH,y,1,dq,1)

! u = y^T*dq*y = a * y

u(:) = a*y(:)

! WorkM = |B|

WorkM(:,:) = abs(H(:,:))

! WorkV = |B|dq

call dGeMV_('N',nH,nH,One,WorkM,nH,dq,1,Zero,WorkV,1)

! b = (dq^T|B|dq) = dq^T * WorkV

b = DDot_(nH,dq,1,WorkV,1)

! u = y^T*dq*y + (dq^T|B|dq)|B|dq = u + b * WorkV

u(:) = u(:)+b*WorkV(:)

! Calculation of f = (y^T*dq)^2 + (dq^T|B|dq)^2 = a^2 + b^2

f = a**2+b**2

! Calculation of v = y - B * dq (quasi-Newton condition)

v(:) = y(:)
call dGeMV_('N',nH,nH,-One,H,nH,dq,1,One,v,1)

! Calculation of equation 23, the new Hessian

! H_k+1 = H_k + 1/f (v*u^T + u*v^T - (a - dq^T*B*dq)/f * u*u^T)

! WorkM = u*u^T

call DGEMM_('N','N',nH,nH,1,One,u,nH,u,1,Zero,WorkM,nH)

! WorkV = dq^T*B

call DGEMM_('N','N',1,nH,nH,One,dq,1,H,nH,Zero,WorkV,1)

! WorkR = (a - dq^T*B*dq)/f = (a - WorkV * dq) / f

WorkR = (a-DDot_(nH,WorkV,1,dq,1))/f

! H_k+1 = H_k + 1/f ( v*u^T + u*v^T - (a - dq^T*B*dq)/f * u*u^T ) =
!         H_k + 1/f ( v*u^T + u*v^T - WorkR * WorkM )

do i=1,nH
  H(:,i) = H(:,i)+One/f*(v(i)*u(:)+u(i)*v(:)-WorkR*WorkM(:,i))
end do

#ifdef _DEBUGPRINT_
! Checking the quasi-Newton condition.

! WorkV = H * dq

call dGeMV_('N',nH,nH,One,H,nH,dq,1,Zero,WorkV,1)

! WorkV = WorkV - y = 0.0

WorkV(:) = WorkV(:)-y(:)
call RecPrt('quasi-Newton',' ',WorkV,1,nH)

write(u6,*) 'good-bye from ts_bfgs'
#endif

call mma_deallocate(WorkM)
call mma_deallocate(WorkV)
call mma_deallocate(v)
call mma_deallocate(u)

return

end subroutine ts_bfgs
