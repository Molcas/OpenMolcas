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

subroutine rotder4(norder,S,X,dXdA,d2XdA2,d3XdA3,d4XdA4)

use Constants, only: Zero, One, Two, Three, Four, Six, Five, Eight, Nine, Ten, Twelve, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: norder
real(kind=wp), intent(in) :: S(3,3), X(3)
real(kind=wp), intent(out) :: dXdA(3,3), d2XdA2(3,3,3), d3XdA3(3,3,3,3), d4XdA4(3,3,3,3,3)
integer(kind=iwp) :: i, i1, i2, i3, ia, ib, ic, id, j, k, l, m, n
real(kind=wp) :: c0, c1, c2, cs, d2AdX2(3,3,3), d2c0dQ2, d2c0dX2(3,3), d2c1dQ2, d2c1dX2(3,3), d2c2dQ2, d2c2dX2(3,3), &
                 d2PdX2(3,3,3,3), d2UdX2(3,3,3,3), d3AdX3(3,3,3,3), d3c0dQ3, d3c0dX3(3,3,3), d3c1dQ3, d3c1dX3(3,3,3), d3c2dQ3, &
                 d3c2dX3(3,3,3), d3PdX3(3,3,3,3,3), d3UdX3(3,3,3,3,3), d4AdX4(3,3,3,3,3), d4c0dQ4, d4c0dX4(3,3,3,3), d4c1dQ4, &
                 d4c1dX4(3,3,3,3), d4c2dQ4, d4c2dX4(3,3,3,3), d4PdX4(3,3,3,3,3,3), d4UdX4(3,3,3,3,3,3), dAdX(3,3), dc0dQ, &
                 dc0dX(3), dc1dQ, dc1dX(3), dc2dQ, dc2dX(3), dPdX(3,3,3), dUdX(3,3,3), Q, S0(3,3), sn, sval(3), tmp1(3), &
                 tmp2(3,3), tmp3(3,3,3), tmp4(3,3,3,3), tmp5(3,3,3,3,3), tmp5A(3,3,3,3,3), U(3,3), umat(3,3), vmat(3,3), &
                 wTmp(100), XN, XX(3,3), XXX(3,3,3), XXXX(3,3,3,3) !, det, detinv, T(3,3)

if (nOrder > 2) call Untested('RotDer4: nOrder > 2')

d2c0dQ2 = Zero ! Dummy initialize
d2c1dQ2 = Zero ! Dummy initialize
d2c2dQ2 = Zero ! Dummy initialize
d3c0dQ3 = Zero ! Dummy initialize
d3c1dQ3 = Zero ! Dummy initialize
d3c2dQ3 = Zero ! Dummy initialize
d4c0dQ4 = Zero ! Dummy initialize
d4c1dQ4 = Zero ! Dummy initialize
d4c2dQ4 = Zero ! Dummy initialize

Q = X(1)**2+X(2)**2+X(3)**2
XN = sqrt(Q)
!                                                                      *
!***********************************************************************
!                                                                      *
if (Q < Half) then
  ! Use a Taylor expansion:
  c0 = One-(Q/Two)*(One-(Q/Twelve)*(One-(Q/30.0_wp)*(One-(Q/56.0_wp))))
  c1 = (One-(Q/Six)*(One-(Q/20.0_wp)*(One-(Q/42.0_wp)*(One-(Q/72.0_wp)))))
  c2 = (One-(Q/Twelve)*(One-(Q/30.0_wp)*(One-(Q/56.0_wp)*(One-(Q/90.0_wp)))))/Two
  dc0dQ = -(One-(Q/Six)*(One-(Q/20.0_wp)*(One-(Q/42.0_wp)*(One-(Q/72.0_wp)))))/Two
  dc1dQ = -(One-(Q/Ten)*(One-(Q/28.0_wp)*(One-(Q/54.0_wp))))/Six
  dc2dQ = -(One-(Q/30.0_wp)*(Two-(Q/56.0_wp)*(Three-(Q/90.0_wp)*Four)))/24.0_wp
  if (norder > 1) then
    d2c0dQ2 = (One-(Q/Ten)*(One-(Q/28.0_wp)*(One-(Q/54.0_wp))))/Twelve
    d2c1dQ2 = (One-(Q/14.0_wp)*(One-(Q/36.0_wp)*(One-(Q/66.0_wp))))/60.0_wp
    d2c2dQ2 = (One-(Q/56.0_wp)*(Three-(Q/90.0_wp)*(Six-(Q/132.0_wp)*Ten)))/360.0_wp
    if (norder > 2) then
      d3c0dQ3 = -(One-(Q/14.0_wp)*(One-(Q/36.0_wp)*(One-(Q/66.0_wp))))/120.0_wp
      d3c1dQ3 = -(One-(Q/18.0_wp)*(One-(Q/44.0_wp)*(One-(Q/78.0_wp))))/840.0_wp
      d3c2dQ3 = -(One-(Q/90.0_wp)*(Four-(Q/132.0_wp)*Ten))/6720.0_wp
      if (norder > 3) then
        d4c0dQ4 = (One-(Q/18.0_wp)*(One-(Q/44.0_wp)*(One-(Q/78.0_wp))))/1680.0_wp
        d4c1dQ4 = (One-(Q/22.0_wp)*(One-(Q/52.0_wp)*(One-(Q/90.0_wp))))/15120.0_wp
        d4c2dQ4 = (One-(Q/132.0_wp)*(Five-(Q/182.0_wp)*15.0_wp))/151200.0_wp
      end if
    end if
  end if
else
  cs = cos(XN)
  sn = sin(XN)
  c0 = cs
  c1 = sn/XN
  c2 = (One-cs)/(XN**2)
  dc0dQ = -sn/(Two*XN)
  dc1dQ = -(sn-XN*cs)/(Two*XN**3)
  dc2dQ = (-Two+Two*cs+XN*sn)/(Two*XN**4)
  if (norder > 1) then
    d2c0dQ2 = (sn-XN*cs)/(Four*XN**3)
    d2c1dQ2 = -(Three*XN*cs-(Three-Q)*sn)/(Four*XN**5)
    d2c2dQ2 = (Eight-(Eight-XN**2)*cs-Five*XN*sn)/(Five*XN**6)
    if (norder > 2) then
      d3c0dQ3 = (Three*XN*cs-(Three-Q)*sn)/(Eight*XN**5)
      d3c1dQ3 = -((15.0_wp-Six*Q)*sn-XN*(15.0_wp-Q)*cs)/(Eight*XN**7)
      d3c2dQ3 = (-48.0_wp+(48.0_wp-Nine*XN**2)*cs+XN*(33.0_wp-XN**2)*sn)/(Eight*XN**8)
      if (norder > 3) then
        d4c0dQ4 = ((15.0_wp-Six*Q)*sn-XN*(15.0_wp-Q)*cs)/(16.0_wp*XN**7)
        d4c1dQ4 = ((105.0_wp-Q*(45.0_wp-Q))*sn-XN*(105.0_wp-Ten*Q)*cs)/(16.0_wp*XN**9)
        d4c2dQ4 = (384.0_wp-(384.0_wp-87.0_wp*XN**2+XN**4)*cs-XN*(279.0_wp-14.0_wp*XN**2)*sn)/(16.0_wp*XN**10)
      end if
    end if
  end if
end if
dc0dX(:) = Two*dc0dQ*X(:)
dc1dX(:) = Two*dc1dQ*X(:)
dc2dX(:) = Two*dc2dQ*X(:)
do j=1,3
  XX(:,j) = X(:)*X(j)
end do
if (norder > 1) then
  d2c0dX2(:,:) = Four*d2c0dQ2*XX(:,:)
  d2c1dX2(:,:) = Four*d2c1dQ2*XX(:,:)
  d2c2dX2(:,:) = Four*d2c2dQ2*XX(:,:)
  do i=1,3
    d2c0dX2(i,i) = d2c0dX2(i,i)+Two*dc0dQ
    d2c1dX2(i,i) = d2c1dX2(i,i)+Two*dc1dQ
    d2c2dX2(i,i) = d2c2dX2(i,i)+Two*dc2dQ
  end do
  if (norder > 2) then
    do k=1,3
      do j=1,3
        XXX(:,j,k) = X(:)*X(j)*X(k)
      end do
    end do
    d3c0dX3(:,:,:) = Eight*d3c0dQ3*XXX(:,:,:)
    d3c1dX3(:,:,:) = Eight*d3c1dQ3*XXX(:,:,:)
    d3c2dX3(:,:,:) = Eight*d3c2dQ3*XXX(:,:,:)
    do i=1,3
      d3c0dX3(i,i,:) = d3c0dX3(i,i,:)+Four*d2c0dQ2*X(:)
      d3c0dX3(i,:,i) = d3c0dX3(i,:,i)+Four*d2c0dQ2*X(:)
      d3c0dX3(:,i,i) = d3c0dX3(:,i,i)+Four*d2c0dQ2*X(:)
      d3c1dX3(i,i,:) = d3c1dX3(i,i,:)+Four*d2c1dQ2*X(:)
      d3c1dX3(i,:,i) = d3c1dX3(i,:,i)+Four*d2c1dQ2*X(:)
      d3c1dX3(:,i,i) = d3c1dX3(:,i,i)+Four*d2c1dQ2*X(:)
      d3c2dX3(i,i,:) = d3c2dX3(i,i,:)+Four*d2c2dQ2*X(:)
      d3c2dX3(i,:,i) = d3c2dX3(i,:,i)+Four*d2c2dQ2*X(:)
      d3c2dX3(:,i,i) = d3c2dX3(:,i,i)+Four*d2c2dQ2*X(:)
    end do
    if (norder > 3) then
      do l=1,3
        do k=1,3
          do j=1,3
            XXXX(:,j,k,l) = X(:)*X(j)*X(k)*X(l)
          end do
        end do
      end do
      d4c0dX4(:,:,:,:) = 16.0_wp*d4c0dQ4*XXXX(:,:,:,:)
      d4c1dX4(:,:,:,:) = 16.0_wp*d4c1dQ4*XXXX(:,:,:,:)
      d4c2dX4(:,:,:,:) = 16.0_wp*d4c2dQ4*XXXX(:,:,:,:)
      do i=1,3
        d4c0dX4(i,i,:,:) = d4c0dX4(i,i,:,:)+Eight*d3c0dQ3*XX(:,:)
        d4c0dX4(i,:,i,:) = d4c0dX4(i,:,i,:)+Eight*d3c0dQ3*XX(:,:)
        d4c0dX4(i,:,:,i) = d4c0dX4(i,:,:,i)+Eight*d3c0dQ3*XX(:,:)
        d4c0dX4(:,i,i,:) = d4c0dX4(:,i,i,:)+Eight*d3c0dQ3*XX(:,:)
        d4c0dX4(:,i,:,i) = d4c0dX4(:,i,:,i)+Eight*d3c0dQ3*XX(:,:)
        d4c0dX4(:,:,i,i) = d4c0dX4(:,:,i,i)+Eight*d3c0dQ3*XX(:,:)
        d4c1dX4(i,i,:,:) = d4c1dX4(i,i,:,:)+Eight*d3c1dQ3*XX(:,:)
        d4c1dX4(i,:,i,:) = d4c1dX4(i,:,i,:)+Eight*d3c1dQ3*XX(:,:)
        d4c1dX4(i,:,:,i) = d4c1dX4(i,:,:,i)+Eight*d3c1dQ3*XX(:,:)
        d4c1dX4(:,i,i,:) = d4c1dX4(:,i,i,:)+Eight*d3c1dQ3*XX(:,:)
        d4c1dX4(:,i,:,i) = d4c1dX4(:,i,:,i)+Eight*d3c1dQ3*XX(:,:)
        d4c1dX4(:,:,i,i) = d4c1dX4(:,:,i,i)+Eight*d3c1dQ3*XX(:,:)
        d4c2dX4(i,i,:,:) = d4c2dX4(i,i,:,:)+Eight*d3c2dQ3*XX(:,:)
        d4c2dX4(i,:,i,:) = d4c2dX4(i,:,i,:)+Eight*d3c2dQ3*XX(:,:)
        d4c2dX4(i,:,:,i) = d4c2dX4(i,:,:,i)+Eight*d3c2dQ3*XX(:,:)
        d4c2dX4(:,i,i,:) = d4c2dX4(:,i,i,:)+Eight*d3c2dQ3*XX(:,:)
        d4c2dX4(:,i,:,i) = d4c2dX4(:,i,:,i)+Eight*d3c2dQ3*XX(:,:)
        d4c2dX4(:,:,i,i) = d4c2dX4(:,:,i,i)+Eight*d3c2dQ3*XX(:,:)
      end do
      do j=1,3
        do i=1,3
          d4c0dX4(i,i,j,j) = d4c0dX4(i,i,j,j)+Four*d2c0dQ2
          d4c0dX4(i,j,i,j) = d4c0dX4(i,j,i,j)+Four*d2c0dQ2
          d4c0dX4(i,j,j,i) = d4c0dX4(i,j,j,i)+Four*d2c0dQ2
          d4c1dX4(i,i,j,j) = d4c1dX4(i,i,j,j)+Four*d2c1dQ2
          d4c1dX4(i,j,i,j) = d4c1dX4(i,j,i,j)+Four*d2c1dQ2
          d4c1dX4(i,j,j,i) = d4c1dX4(i,j,j,i)+Four*d2c1dQ2
          d4c2dX4(i,i,j,j) = d4c2dX4(i,i,j,j)+Four*d2c2dQ2
          d4c2dX4(i,j,i,j) = d4c2dX4(i,j,i,j)+Four*d2c2dQ2
          d4c2dX4(i,j,j,i) = d4c2dX4(i,j,j,i)+Four*d2c2dQ2
        end do
      end do
    end if
  end if
end if
! The unitary matrix U=exp(XMat), and its derivatives.
! Here, X is the matrix with elements XMat(i,j)=eps(i,k,j)*X(k)
! First term, cos(X)*delta(i,j):
U(:,:) = Zero
do i=1,3
  U(i,i) = c0
end do
dUdX(:,:,:) = Zero
do i=1,3
  dUdX(i,i,:) = dc0dX(:)
end do

if (norder > 1) then
  d2UdX2(:,:,:,:) = Zero
  do i=1,3
    d2UdX2(i,i,:,:) = d2c0dX2(:,:)
  end do
  if (norder > 2) then
    d3UdX3(:,:,:,:,:) = Zero
    do i=1,3
      d3UdX3(i,i,:,:,:) = d3c0dX3(:,:,:)
    end do
    if (norder > 3) then
      d4UdX4(:,:,:,:,:,:) = Zero
      do i=1,3
        d4UdX4(i,i,:,:,:,:) = d4c0dX4(:,:,:,:)
      end do
    end if
  end if
end if
! Second term, (sin(X)/X)*eps(i,k,j)*X(k):
do i1=1,3
  i2 = 1+mod(i1,3)
  i3 = 1+mod(i2,3)
  U(i1,i3) = U(i1,i3)+c1*X(i2)
  U(i3,i1) = U(i3,i1)-c1*X(i2)
  dUdX(i1,i3,:) = dUdX(i1,i3,:)+dc1dX(:)*X(i2)
  dUdX(i3,i1,:) = dUdX(i3,i1,:)-dc1dX(:)*X(i2)
  dUdX(i1,i3,i2) = dUdX(i1,i3,i2)+c1
  dUdX(i3,i1,i2) = dUdX(i3,i1,i2)-c1
  if (norder > 1) then
    d2UdX2(i1,i3,:,:) = d2UdX2(i1,i3,:,:)+d2c1dX2(:,:)*X(i2)
    d2UdX2(i3,i1,:,:) = d2UdX2(i3,i1,:,:)-d2c1dX2(:,:)*X(i2)
    d2UdX2(i1,i3,i2,:) = d2UdX2(i1,i3,i2,:)+dc1dX(:)
    d2UdX2(i3,i1,i2,:) = d2UdX2(i3,i1,i2,:)-dc1dX(:)
    d2UdX2(i1,i3,:,i2) = d2UdX2(i1,i3,:,i2)+dc1dX(:)
    d2UdX2(i3,i1,:,i2) = d2UdX2(i3,i1,:,i2)-dc1dX(:)
    if (norder > 2) then
      d3UdX3(i1,i3,:,:,:) = d3UdX3(i1,i3,:,:,:)+d3c1dX3(:,:,:)*X(i2)
      d3UdX3(i3,i1,:,:,:) = d3UdX3(i3,i1,:,:,:)-d3c1dX3(:,:,:)*X(i2)
      d3UdX3(i1,i3,i2,:,:) = d3UdX3(i1,i3,i2,:,:)+d2c1dX2(:,:)
      d3UdX3(i3,i1,i2,:,:) = d3UdX3(i3,i1,i2,:,:)-d2c1dX2(:,:)
      d3UdX3(i1,i3,:,i2,:) = d3UdX3(i1,i3,:,i2,:)+d2c1dX2(:,:)
      d3UdX3(i3,i1,:,i2,:) = d3UdX3(i3,i1,:,i2,:)-d2c1dX2(:,:)
      d3UdX3(i1,i3,:,:,i2) = d3UdX3(i1,i3,:,:,i2)+d2c1dX2(:,:)
      d3UdX3(i3,i1,:,:,i2) = d3UdX3(i3,i1,:,:,i2)-d2c1dX2(:,:)
      if (norder > 3) then
        d4UdX4(i1,i3,:,:,:,:) = d4UdX4(i1,i3,:,:,:,:)+d4c1dX4(:,:,:,:)*X(i2)
        d4UdX4(i3,i1,:,:,:,:) = d4UdX4(i3,i1,:,:,:,:)-d4c1dX4(:,:,:,:)*X(i2)
        d4UdX4(i1,i3,i2,:,:,:) = d4UdX4(i1,i3,i2,:,:,:)+d3c1dX3(:,:,:)
        d4UdX4(i3,i1,i2,:,:,:) = d4UdX4(i3,i1,i2,:,:,:)-d3c1dX3(:,:,:)
        d4UdX4(i1,i3,:,i2,:,:) = d4UdX4(i1,i3,:,i2,:,:)+d3c1dX3(:,:,:)
        d4UdX4(i3,i1,:,i2,:,:) = d4UdX4(i3,i1,:,i2,:,:)-d3c1dX3(:,:,:)
        d4UdX4(i1,i3,:,:,i2,:) = d4UdX4(i1,i3,:,:,i2,:)+d3c1dX3(:,:,:)
        d4UdX4(i3,i1,:,:,i2,:) = d4UdX4(i3,i1,:,:,i2,:)-d3c1dX3(:,:,:)
        d4UdX4(i1,i3,:,:,:,i2) = d4UdX4(i1,i3,:,:,:,i2)+d3c1dX3(:,:,:)
        d4UdX4(i3,i1,:,:,:,i2) = d4UdX4(i3,i1,:,:,:,i2)-d3c1dX3(:,:,:)
      end if
    end if
  end if
end do
! Third term,((1-cos(X))/X**2)*X(i)*X(j):
U(:,:) = U(:,:)+c2*XX(:,:)
do k=1,3
  dUdX(:,:,k) = dUdX(:,:,k)+dc2dX(k)*XX(:,:)
end do
do i=1,3
  dUdX(i,:,i) = dUdX(i,:,i)+c2*X(:)
end do
do j=1,3
  dUdX(:,j,j) = dUdX(:,j,j)+c2*X(:)
end do
if (norder > 1) then
  do l=1,3
    do k=1,3
      d2UdX2(:,:,k,l) = d2UdX2(:,:,k,l)+d2c2dX2(k,l)*XX(:,:)
    end do
  end do
  do k=1,3
    do i=1,3
      d2UdX2(i,:,i,k) = d2UdX2(i,:,i,k)+dc2dX(k)*X(:)
      d2UdX2(i,:,k,i) = d2UdX2(i,:,k,i)+dc2dX(k)*X(:)
    end do
    do j=1,3
      d2UdX2(:,j,j,k) = d2UdX2(:,j,j,k)+dc2dX(k)*X(:)
      d2UdX2(:,j,k,j) = d2UdX2(:,j,k,j)+dc2dX(k)*X(:)
    end do
  end do
  do j=1,3
    do i=1,3
      d2UdX2(i,j,i,j) = d2UdX2(i,j,i,j)+c2
      d2UdX2(i,j,j,i) = d2UdX2(i,j,j,i)+c2
    end do
  end do
  if (norder > 2) then
    do m=1,3
      do l=1,3
        do k=1,3
          d3UdX3(:,:,k,l,m) = d3UdX3(:,:,k,l,m)+d3c2dX3(k,l,m)*XX(:,:)
        end do
      end do
    end do
    do l=1,3
      do k=1,3
        do i=1,3
          d3UdX3(i,:,i,k,l) = d3UdX3(i,:,i,k,l)+d2c2dX2(k,l)*X(:)
          d3UdX3(i,:,k,i,l) = d3UdX3(i,:,k,i,l)+d2c2dX2(k,l)*X(:)
          d3UdX3(i,:,k,l,i) = d3UdX3(i,:,k,l,i)+d2c2dX2(k,l)*X(:)
        end do
        do j=1,3
          d3UdX3(:,j,j,k,l) = d3UdX3(:,j,j,k,l)+d2c2dX2(k,l)*X(:)
          d3UdX3(:,j,k,j,l) = d3UdX3(:,j,k,j,l)+d2c2dX2(k,l)*X(:)
          d3UdX3(:,j,k,l,j) = d3UdX3(:,j,k,l,j)+d2c2dX2(k,l)*X(:)
        end do
      end do
    end do
    do j=1,3
      do i=1,3
        d3UdX3(i,j,i,j,:) = d3UdX3(i,j,i,j,:)+dc2dX(:)
        d3UdX3(i,j,j,i,:) = d3UdX3(i,j,j,i,:)+dc2dX(:)
        d3UdX3(i,j,i,:,j) = d3UdX3(i,j,i,:,j)+dc2dX(:)
        d3UdX3(i,j,j,:,i) = d3UdX3(i,j,j,:,i)+dc2dX(:)
        d3UdX3(i,j,:,i,j) = d3UdX3(i,j,:,i,j)+dc2dX(:)
        d3UdX3(i,j,:,j,i) = d3UdX3(i,j,:,j,i)+dc2dX(:)
      end do
    end do
    if (norder > 3) then
      do n=1,3
        do m=1,3
          do l=1,3
            do k=1,3
              d4UdX4(:,:,k,l,m,n) = d4UdX4(:,:,k,l,m,n)+d4c2dX4(k,l,m,n)*XX(:,:)
            end do
          end do
        end do
      end do
      do m=1,3
        do l=1,3
          do k=1,3
            do i=1,3
              d4UdX4(i,:,i,k,l,m) = d4UdX4(i,:,i,k,l,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(i,:,k,i,l,m) = d4UdX4(i,:,k,i,l,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(i,:,k,l,i,m) = d4UdX4(i,:,k,l,i,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(i,:,k,l,m,i) = d4UdX4(i,:,k,l,m,i)+d3c2dX3(k,l,m)*X(:)
            end do
            do j=1,3
              d4UdX4(:,j,j,k,l,m) = d4UdX4(:,j,j,k,l,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(:,j,k,j,l,m) = d4UdX4(:,j,k,j,l,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(:,j,k,l,j,m) = d4UdX4(:,j,k,l,j,m)+d3c2dX3(k,l,m)*X(:)
              d4UdX4(:,j,k,l,m,j) = d4UdX4(:,j,k,l,m,j)+d3c2dX3(k,l,m)*X(:)
            end do
          end do
        end do
      end do
      do j=1,3
        do i=1,3
          d4UdX4(i,j,i,j,:,:) = d4UdX4(i,j,i,j,:,:)+d2c2dX2(:,:)
          d4UdX4(i,j,i,:,j,:) = d4UdX4(i,j,i,:,j,:)+d2c2dX2(:,:)
          d4UdX4(i,j,i,:,:,j) = d4UdX4(i,j,i,:,:,j)+d2c2dX2(:,:)
          d4UdX4(i,j,:,i,j,:) = d4UdX4(i,j,:,i,j,:)+d2c2dX2(:,:)
          d4UdX4(i,j,:,i,:,j) = d4UdX4(i,j,:,i,:,j)+d2c2dX2(:,:)
          d4UdX4(i,j,:,:,i,j) = d4UdX4(i,j,:,:,i,j)+d2c2dX2(:,:)
          d4UdX4(i,j,j,i,:,:) = d4UdX4(i,j,j,i,:,:)+d2c2dX2(:,:)
          d4UdX4(i,j,j,:,i,:) = d4UdX4(i,j,j,:,i,:)+d2c2dX2(:,:)
          d4UdX4(i,j,j,:,:,i) = d4UdX4(i,j,j,:,:,i)+d2c2dX2(:,:)
          d4UdX4(i,j,:,j,i,:) = d4UdX4(i,j,:,j,i,:)+d2c2dX2(:,:)
          d4UdX4(i,j,:,j,:,i) = d4UdX4(i,j,:,j,:,i)+d2c2dX2(:,:)
          d4UdX4(i,j,:,:,j,i) = d4UdX4(i,j,:,:,j,i)+d2c2dX2(:,:)
        end do
      end do
    end if
  end if
end if
! The matrix S (which should be symmetrical) is the product
! S=S0*U=S0*exp(X) with X given as input arguments to this Subroutine.
! Need to compute S0=S*U(transpose):
S0(:,:) = Zero
do i=1,3
  do j=1,3
    do k=1,3
      S0(i,j) = S0(i,j)+S(i,k)*U(j,k)
    end do
  end do
end do
! The matrix product P=S0*U, where now U is used as an expansion in
! the variables X, and its derivatives:
do i=1,3
  do j=1,3
    !P(i,j) = sum(S0(i,:)*U(:,j))
    do k=1,3
      dPdX(i,j,k) = sum(S0(i,:)*dUdX(:,j,k))
    end do
  end do
end do
if (norder > 1) then
  do l=1,3
    do k=1,3
      do j=1,3
        do i=1,3
          d2PdX2(i,j,k,l) = sum(S0(i,:)*d2UdX2(:,j,k,l))
        end do
      end do
    end do
  end do
  if (norder > 2) then
    do m=1,3
      do l=1,3
        do k=1,3
          do j=1,3
            do i=1,3
              d3PdX3(i,j,k,l,m) = sum(S0(i,:)*d3UdX3(:,j,k,l,m))
            end do
          end do
        end do
      end do
    end do
    if (norder > 3) then
      do n=1,3
        do m=1,3
          do l=1,3
            do k=1,3
              do j=1,3
                do i=1,3
                  d4PdX4(i,j,k,l,m,n) = sum(S0(i,:)*d4UdX4(:,j,k,l,m,n))
                end do
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if
! The vector A is the dual of the antisymmetric part of P:
!       A(1) = P(3,2)-P(2,3)
!       A(2) = P(1,3)-P(3,1)
!       A(3) = P(2,1)-P(1,2)
dAdX(1,:) = dPdX(3,2,:)-dPdX(2,3,:)
dAdX(2,:) = dPdX(1,3,:)-dPdX(3,1,:)
dAdX(3,:) = dPdX(2,1,:)-dPdX(1,2,:)
if (norder > 1) then
  d2AdX2(1,:,:) = d2PdX2(3,2,:,:)-d2PdX2(2,3,:,:)
  d2AdX2(2,:,:) = d2PdX2(1,3,:,:)-d2PdX2(3,1,:,:)
  d2AdX2(3,:,:) = d2PdX2(2,1,:,:)-d2PdX2(1,2,:,:)
  if (norder > 2) then
    d3AdX3(1,:,:,:) = d3PdX3(3,2,:,:,:)-d3PdX3(2,3,:,:,:)
    d3AdX3(2,:,:,:) = d3PdX3(1,3,:,:,:)-d3PdX3(3,1,:,:,:)
    d3AdX3(3,:,:,:) = d3PdX3(2,1,:,:,:)-d3PdX3(1,2,:,:,:)
    if (norder > 3) then
      d4AdX4(1,:,:,:,:) = d4PdX4(3,2,:,:,:,:)-d4PdX4(2,3,:,:,:,:)
      d4AdX4(2,:,:,:,:) = d4PdX4(1,3,:,:,:,:)-d4PdX4(3,1,:,:,:,:)
      d4AdX4(3,:,:,:,:) = d4PdX4(2,1,:,:,:,:)-d4PdX4(1,2,:,:,:,:)
    end if
  end if
end if
! Finally, we have obtained derivatives of A w.r.t X.
! So now we can obtain derivatives of the inverse mapping.
! This will require the inverse T of the matrix dAdX:
!     T(1,1) = dAdX(2,2)*dAdX(3,3)-dAdX(3,2)*dAdX(2,3)
!     T(2,1) = dAdX(3,2)*dAdX(1,3)-dAdX(1,2)*dAdX(3,3)
!     T(3,1) = dAdX(1,2)*dAdX(2,3)-dAdX(2,2)*dAdX(1,3)
!     T(1,2) = dAdX(2,3)*dAdX(3,1)-dAdX(3,3)*dAdX(2,1)
!     T(2,2) = dAdX(3,3)*dAdX(1,1)-dAdX(1,3)*dAdX(3,1)
!     T(3,2) = dAdX(1,3)*dAdX(2,1)-dAdX(2,3)*dAdX(1,1)
!     T(1,3) = dAdX(2,1)*dAdX(3,2)-dAdX(3,1)*dAdX(2,2)
!     T(2,3) = dAdX(3,1)*dAdX(1,2)-dAdX(1,1)*dAdX(3,2)
!     T(3,3) = dAdX(1,1)*dAdX(2,2)-dAdX(2,1)*dAdX(1,2)
!     det = dAdX(1,1)*T(1,1)+dAdX(2,1)*T(2,1)+dAdX(3,1)*T(3,1)
!     detInv = One/det
! First derivatives dXdA(j,i):
!     do i=1,3
!       do ia=1,3
!         dXdA(ia,i) = DetInv*T(i,ia)
!       end do
!     end do
! Use the Moore-Penrose pseudoinverse instead
call dgesvd_('A','A',3,3,dAdX,3,sval,umat,3,vmat,3,wTmp,100,i)
do i=1,3
  if (abs(sval(i)) > 1.0e-12_wp) then
    umat(:,i) = umat(:,i)/sval(i)
  else
    umat(:,i) = Zero
  end if
end do
call dgemm_('T','T',3,3,3,One,vmat,3,umat,3,Zero,dXdA,3)
! Second derivatives d2XdA(ic,j,k)
if (norder > 1) then
  do k=1,3
    do i=1,3
      do ia=1,3
        tmp1(ia) = sum(d2AdX2(i,ia,:)*dXdA(:,k))
      end do
      do j=1,3
        tmp3(i,j,k) = sum(tmp1(:)*dXdA(:,j))
      end do
    end do
  end do
  do k=1,3
    do j=1,3
      do ic=1,3
        d2XdA2(ic,j,k) = -sum(dXdA(ic,:)*tmp3(:,j,k))
      end do
    end do
  end do
  ! Third derivatives d3XdA3(id,j,k,l)
  if (norder > 2) then
    do l=1,3
      do i=1,3
        do ia=1,3
          do ib=1,3
            tmp1(ib) = sum(d3AdX3(i,ia,ib,:)*dXdA(:,l))
          end do
          do k=1,3
            tmp2(ia,k) = sum(tmp1(:)*dXdA(:,k))
          end do
        end do
        do k=1,3
          do j=1,3
            tmp4(i,j,k,l) = sum(tmp2(:,k)*dXdA(:,j))
          end do
        end do
      end do
    end do
    do l=1,3
      do k=1,3
        do j=1,3
          do id=1,3
            d3XdA3(id,j,k,l) = -sum(dXdA(id,:)*tmp4(:,j,k,l))
          end do
        end do
      end do
    end do
    do j=1,3
      do i=1,3
        do ib=1,3
          tmp1(ib) = sum(d2AdX2(i,:,ib)*dXdA(:,j))
        end do
        do l=1,3
          do k=1,3
            tmp4(i,j,k,l) = sum(tmp1(:)*d2XdA2(:,k,l))
          end do
        end do
      end do
    end do
    do l=1,3
      do k=1,3
        do j=1,3
          do ic=1,3
            d3XdA3(ic,j,k,l) = d3XdA3(ic,j,k,l)-sum(DXDA(ic,:)*(tmp4(:,j,k,l)+tmp4(:,k,l,j)+tmp4(:,l,j,k)))
          end do
        end do
      end do
    end do
    ! Fourth derivatives d4XdA4(id,j,k,l,m)
    if (norder > 3) then
      do i=1,3
        do ia=1,3
          do ib=1,3
            do m=1,3
              do ic=1,3
                tmp2(ic,m) = sum(d4AdX4(i,ia,ib,ic,:)*dXdA(:,m))
              end do
            end do
            do m=1,3
              do l=1,3
                tmp3(ib,l,m) = sum(tmp2(:,m)*dXdA(:,l))
              end do
            end do
          end do
          do m=1,3
            do l=1,3
              do k=1,3
                tmp4(ia,k,l,m) = sum(tmp3(:,l,m)*dXdA(:,k))
              end do
            end do
          end do
        end do
        do m=1,3
          do l=1,3
            do k=1,3
              do j=1,3
                tmp5(i,j,k,l,m) = sum(tmp4(:,k,l,m)*dXdA(:,j))
              end do
            end do
          end do
        end do
        do l=1,3
          do m=1,3
            do ia=1,3
              do ib=1,3
                tmp1(ib) = sum(d3AdX3(i,ia,ib,:)*d2XdA2(:,l,m))
              end do
              do k=1,3
                tmp2(ia,k) = sum(tmp1(:)*dXdA(:,k))
              end do
            end do
            do k=1,3
              do j=1,3
                tmp5A(i,j,k,l,m) = sum(tmp2(:,k)*dXdA(:,j))
              end do
            end do
          end do
        end do
        do m=1,3
          do l=1,3
            do k=1,3
              do j=1,3
                tmp5(i,j,k,l,m) = tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+tmp5A(i,j,l,k,m)+tmp5A(i,k,l,j,m)+tmp5A(i,j,m,k,l)+ &
                                  tmp5A(i,k,m,j,l)+tmp5A(i,l,m,j,k)
              end do
            end do
          end do
        end do
        do j=1,3
          do ib=1,3
            tmp2(j,ib) = sum(d2AdX2(i,:,ib)*dXdA(:,j))
          end do
          do m=1,3
            do l=1,3
              do k=1,3
                tmp5A(i,j,k,l,m) = sum(tmp2(j,:)*d3XdA3(:,k,l,m))
              end do
            end do
          end do
        end do
        do m=1,3
          do l=1,3
            do k=1,3
              do j=1,3
                tmp5(i,j,k,l,m) = tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+tmp5A(i,k,l,m,j)+tmp5A(i,l,m,j,k)+tmp5A(i,m,j,k,l)
              end do
            end do
          end do
        end do
        do j=1,3
          do m=1,3
            do ib=1,3
              tmp1(ib) = sum(d2AdX2(i,:,ib)*d2XdA2(:,j,m))
            end do
            do l=1,3
              do k=1,3
                tmp5A(i,j,k,l,m) = sum(tmp1(:)*d2XdA2(:,k,l))
              end do
            end do
          end do
        end do
        do m=1,3
          do l=1,3
            do k=1,3
              do j=1,3
                tmp5(i,j,k,l,m) = tmp5(i,j,k,l,m)+tmp5A(i,j,k,l,m)+tmp5A(i,j,k,m,l)+tmp5A(i,j,m,l,k)
              end do
            end do
          end do
        end do
      end do
      do m=1,3
        do l=1,3
          do k=1,3
            do j=1,3
              do id=1,3
                d4XdA4(id,j,k,l,m) = -sum(dXdA(id,:)*tmp5(:,j,k,l,m))
              end do
            end do
          end do
        end do
      end do
    end if
  end if
end if

return

end subroutine rotder4
