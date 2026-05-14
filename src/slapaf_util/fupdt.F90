!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine FUpdt(nInter,FEq,gi_2,gi_1,g_ref,qi_2,qi_1,q_ref,u,v,w)
!***********************************************************************
!                                                                      *
! Object: to do a rank-three update on the anharmonic constants.       *
!                                                                      *
!         Observe that SlapAf stores the forces rather than the        *
!         gradients.                                                   *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nInter
real(kind=wp), intent(inout) :: FEq(nInter,nInter,nInter)
real(kind=wp), intent(in) :: gi_2(nInter), gi_1(nInter), g_ref(nInter), qi_2(nInter), qi_1(nInter), q_ref(nInter)
real(kind=wp), intent(out) :: u(nInter), v(nInter), w(nInter)
integer(kind=iwp) :: i, k, l
real(kind=wp) :: lambda, RHS, rLHS, ux, uy, vx, vy, wx, wy
real(kind=wp), external :: DDot_

!call RecPrt(' FEq',' ',FEq,nInter**2,nInter)
!call RecPrt('  gi-2',' ', gi_2,1,nInter)
!call RecPrt('  gi-1',' ', gi_1,1,nInter)
!call RecPrt('  g_ref',' ', g_ref,1,nInter)
!call RecPrt('  qi-2',' ', qi_2,1,nInter)
!call RecPrt('  qi-1',' ', qi_1,1,nInter)
!call RecPrt('  q_ref',' ', q_ref,1,nInter)

u(:) = -(gi_1(:)-g_ref(:))
v(:) = -(gi_2(:)-g_ref(:))
rLHS = DDot_(nInter,qi_2,1,u,1)-DDot_(nInter,q_ref,1,u,1)-DDot_(nInter,qi_1,1,v,1)+DDot_(nInter,q_ref,1,v,1)
write(u6,*) 'FUpdt: LHS=',rLHS
RHS = Zero
do i=1,nInter
  do k=1,nInter
    do l=1,nInter
      RHS = RHS+FEq(i,k,l)*(qi_1(i)-q_ref(i))*(qi_2(k)-q_ref(k))*(qi_2(l)-qi_1(l))
    end do
  end do
end do
RHS = Half*RHS
write(u6,*) 'FUpdt: RHS=',RHS
lambda = rLHS-RHS
write(u6,*) ' FUpdt: lambda=',lambda
w(:) = v(:)-u(:)
call RecPrt('u',' ',u,1,nInter)
call RecPrt('v',' ',v,1,nInter)
call RecPrt('w',' ',w,1,nInter)

ux = DDot_(nInter,u,1,qi_1,1)-DDot_(nInter,u,1,q_ref,1)
uy = DDot_(nInter,u,1,qi_2,1)-DDot_(nInter,u,1,q_ref,1)
vx = DDot_(nInter,v,1,qi_1,1)-DDot_(nInter,v,1,q_ref,1)
vy = DDot_(nInter,v,1,qi_2,1)-DDot_(nInter,v,1,q_ref,1)
wx = DDot_(nInter,w,1,qi_1,1)-DDot_(nInter,w,1,q_ref,1)
wy = DDot_(nInter,w,1,qi_2,1)-DDot_(nInter,w,1,q_ref,1)
lambda = Two*lambda/(ux*vy*(wy-wx)+vx*wy*(uy-ux)+wx*uy*(vy-vx))
write(u6,*) ' FUpdt: lambda=',lambda
do l=1,nInter
  do k=1,nInter
    FEq(:,k,l) = FEq(:,k,l)+lambda*(u(:)*v(k)*w(l)+v(:)*w(k)*u(l)+w(:)*u(k)*v(l))
  end do
end do

! Check on the third order condition

u(:) = -(gi_1(:)-g_ref(:))
v(:) = -(gi_2(:)-g_ref(:))
rLHS = DDot_(nInter,qi_2,1,u,1)-DDot_(nInter,q_ref,1,u,1)-DDot_(nInter,qi_1,1,v,1)+DDot_(nInter,q_ref,1,v,1)
write(u6,*) 'FUpdt: LHS(qNR)=',rLHS
RHS = Zero
do i=1,nInter
  do k=1,nInter
    do l=1,nInter
      RHS = RHS+FEq(i,k,l)*(qi_1(i)-q_ref(i))*(qi_2(k)-q_ref(k))*(qi_2(l)-qi_1(l))
    end do
  end do
end do
RHS = Half*RHS
write(u6,*) 'FUpdt: RHS(qNR)=',RHS

return

end subroutine FUpdt
