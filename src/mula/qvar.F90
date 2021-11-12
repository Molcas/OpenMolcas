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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine qvar_to_var(var,x,grad,Hess,D3,D4,ref,qref,trfName,alpha,max_term,ndata,nvar)
!  Purpose:
!    Tranform coordinates, gradient, Hessian, third derivatives and
!    fourth derivates back to the coordinates originally specified.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use Constants, only: Zero, One, Three, Four, Six
use Definitions, only: u6

implicit real*8(a-h,o-z)
real*8 var(ndata,nvar)
real*8 x(nvar)
real*8 q(nvar)
real*8 t(nvar), u(nvar), v(nvar), s(nvar)
real*8 alpha(nvar)
real*8 grad(nvar), tempgrad(nvar)
real*8 Hess(nvar,nvar)
real*8 Temp(nvar,nvar)
real*8 D3(nvar,nvar,nvar)
real*8 D3trans(nvar,nvar,nvar)
real*8 D4(nvar,nvar,nvar,nvar)
real*8 D4trans(nvar,nvar,nvar,nvar)
real*8 ref(nvar), qref(nvar)
!integer nOrd(4)
character*80 trfName(nvar)
character*32 trfCode
logical ijEq, jkEq, klEq

x(:) = x+qref
do ivar=1,nvar
  trfcode = trfName(ivar)(1:32)
  ix = index(trfcode,'AS IT IS')
  ia = index(trfcode,'-AVG')
  ie = index(trfcode,'EXP')
  ic = index(trfcode,'COS')
  is = index(trfcode,'SIN')
  if (ic > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)-ref(ivar)
    end if
    q(ivar) = acos(-x(ivar))
    t(ivar) = sin(q(ivar))
    u(ivar) = cos(q(ivar))
    v(ivar) = -sin(q(ivar))
    s(ivar) = -cos(q(ivar))
  else if (is > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)-ref(ivar)
    end if
    q(ivar) = asin(-x(ivar))
    t(ivar) = -cos(q(ivar))
    u(ivar) = sin(q(ivar))
    v(ivar) = cos(q(ivar))
    s(ivar) = -sin(q(ivar))
  else if (ie > 0) then
    q(ivar) = -log(One-x(ivar))/alpha(ivar)
    t(ivar) = alpha(ivar)*(One-x(ivar))
    u(ivar) = -alpha(ivar)**2*(One-x(ivar))
    v(ivar) = alpha(ivar)**3*(One-x(ivar))
    s(ivar) = -alpha(ivar)**4*(One-x(ivar))
    if (ia > 0) then
      q(ivar) = q(ivar)+ref(ivar)
    end if
  else if (ix > 0) then
    if (ia > 0) then
      x(ivar) = x(ivar)+ref(ivar)
    end if
    q(ivar) = x(ivar)
    t(ivar) = One
    u(ivar) = Zero
    v(ivar) = Zero
    s(ivar) = Zero
  else
    write(u6,*) ' TRFCODE ERROR.'
    call abend()
  end if
end do

! Transform gradient.
do i=1,nvar
  tempgrad(i) = t(i)*grad(i)
end do

! Transform Hessian.
do i=1,nvar
  do j=i,nvar
    if (i == j) then
      Temp(i,i) = t(i)**2*Hess(i,i)+u(i)*grad(i)
    else
      Temp(i,j) = t(i)*t(j)*Hess(i,j)
    end if
  end do
end do
do i=1,nvar
  do j=i,nvar
    Temp(j,i) = Temp(i,j)
  end do
end do

! Transform third derivatives.
do i=1,nvar
  do j=i,nvar
    ijEq = (i == j)
    do k=j,nvar
      jkEq = (j == k)
      if (ijEq .and. jkEq) then
        D3trans(i,i,i) = t(i)*t(i)*t(i)*D3(i,i,i)+Three*u(i)*t(i)*Hess(i,i)+v(i)*grad(i)
      else if (ijEq .and. (.not. jkEq)) then
        D3trans(i,i,k) = t(i)*t(i)*t(k)*D3(i,i,k)+u(i)*t(k)*Hess(i,k)
      else if ((.not. ijEq) .and. jkEq) then
        D3trans(i,j,j) = t(i)*t(j)*t(j)*D3(i,j,j)+t(i)*u(j)*Hess(i,j)
      else if ((.not. ijEq) .and. (.not. jkEq)) then
        D3trans(i,j,k) = t(i)*t(j)*t(k)*D3(i,j,k)
      end if
    end do
  end do
end do
do i=1,nvar
  do j=i,nvar
    do k=j,nvar
      D3trans(i,k,j) = D3trans(i,j,k)
      D3trans(j,i,k) = D3trans(i,j,k)
      D3trans(j,k,i) = D3trans(i,j,k)
      D3trans(k,i,j) = D3trans(i,j,k)
      D3trans(k,j,i) = D3trans(i,j,k)
    end do
  end do
end do

! Transform fourth derivatives.
do i=1,nvar
  do j=i,nvar
    ijEq = (i == j)
    do k=j,nvar
      jkEq = (j == k)
      do l=k,nvar
        klEq = (k == l)
        if (ijEq .and. jkEq .and. klEq) then
          D4trans(i,i,i,i) = t(i)**4*D4(i,i,i,i)+Six*u(i)*t(i)**2*D3(i,i,i)+Four*v(i)*t(i)*Hess(i,i)+Three*u(i)*u(i)*Hess(i,i)+ &
                             s(i)*grad(i)
        else if (ijEq .and. jkEq .and. (.not. klEq)) then
          D4trans(i,i,i,l) = (t(i)*t(i)*t(i)*D4(i,i,i,l)+Three*u(i)*t(i)*D3(i,i,l)+v(i)*Hess(i,l))*t(l)
        else if ((.not. ijEq) .and. jkEq .and. klEq) then
          D4trans(i,j,j,j) = t(i)*(t(j)*t(j)*t(j)*D4(i,j,j,j)+Three*u(j)*t(j)*D3(i,j,j)+v(j)*Hess(i,j))
        else if ((.not. ijEq) .and. jkEq .and. (.not. klEq)) then
          D4trans(i,j,j,l) = t(i)*(t(j)*t(j)*D4(i,j,j,l)+u(j)*D3(i,j,l))*t(l)
        else if (ijEq .and. (.not. jkEq) .and. (.not. klEq)) then
          D4trans(i,i,k,l) = (t(i)*t(i)*D4(i,i,k,l)+u(i)*D3(i,k,l))*t(k)*t(l)
        else if ((.not. ijEq) .and. (.not. jkEq) .and. klEq) then
          D4trans(i,j,k,k) = (t(k)*t(k)*D4(i,j,k,k)+u(k)*D3(i,j,k))*t(i)*t(j)
        else if (ijEq .and. klEq) then
          D4trans(i,i,k,k) = t(i)*t(i)*t(k)*t(k)*D4(i,i,k,k)+t(i)*t(i)*u(k)*D3(i,i,k)+u(i)*t(k)*t(k)*D3(i,k,k)+u(i)*u(k)*Hess(i,k)
        else if ((.not. ijEq) .and. (.not. jkEq) .and. (.not. klEq)) then
          D4trans(i,j,k,l) = t(i)*t(j)*t(k)*t(l)*D4(i,j,k,l)
        end if
      end do
    end do
  end do
end do
do i=1,nvar
  do j=i,nvar
    do k=j,nvar
      do l=k,nvar
        D4trans(i,k,j,l) = D4trans(i,j,k,l)
        D4trans(i,j,l,k) = D4trans(i,j,k,l)
        D4trans(i,l,k,j) = D4trans(i,j,k,l)
        D4trans(i,l,j,k) = D4trans(i,j,k,l)
        D4trans(i,k,l,j) = D4trans(i,j,k,l)

        D4trans(j,i,k,l) = D4trans(i,j,k,l)
        D4trans(j,k,i,l) = D4trans(i,j,k,l)
        D4trans(j,i,l,k) = D4trans(i,j,k,l)
        D4trans(j,l,k,i) = D4trans(i,j,k,l)
        D4trans(j,l,i,k) = D4trans(i,j,k,l)
        D4trans(j,k,l,i) = D4trans(i,j,k,l)

        D4trans(k,i,j,l) = D4trans(i,j,k,l)
        D4trans(k,j,i,l) = D4trans(i,j,k,l)
        D4trans(k,i,l,j) = D4trans(i,j,k,l)
        D4trans(k,l,j,i) = D4trans(i,j,k,l)
        D4trans(k,l,i,j) = D4trans(i,j,k,l)
        D4trans(k,j,l,i) = D4trans(i,j,k,l)

        D4trans(l,i,j,k) = D4trans(i,j,k,l)
        D4trans(l,j,i,k) = D4trans(i,j,k,l)
        D4trans(l,k,j,i) = D4trans(i,j,k,l)
        D4trans(l,i,k,j) = D4trans(i,j,k,l)
        D4trans(l,k,i,j) = D4trans(i,j,k,l)
        D4trans(l,j,k,i) = D4trans(i,j,k,l)
      end do
    end do
  end do
end do

! Assign transformed values.
x(:) = q
grad(:) = tempgrad

Hess(:,:) = Temp

D3(:,:,:) = D3trans

D4(:,:,:,:) = D4trans

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(var)
  call Unused_integer(max_term)
end if

end subroutine qvar_to_var
