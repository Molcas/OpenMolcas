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

subroutine hotred(nx,n,a,d,e,z)

implicit real*8(a-h,o-z)
dimension a(nx*(nx+1)/2), z(nx,nx), d(nx), e(nx)

if (n <= 2) then
  select case (n)
    case default ! (1)
      d(1) = a(1)
      e(1) = 0.0d0
      z(1,1) = 1.0d0
    case (2)
      d(1) = a(1)
      d(2) = a(3)
      e(1) = 0.0d0
      e(2) = a(2)
      z(1,1) = 1.0d0
      z(2,2) = 1.0d0
      z(1,2) = 0.0d0
      z(2,1) = 0.0d0
  end select
  return
end if

ij = 0
do i=1,n
  do j=1,i
    ij = ij+1
    z(i,j) = a(ij)
  end do
end do
do ip=2,n
  i = n-ip+2
  l = i-2
  f = z(i,i-1)
  g = 0.0d0
  if (l /= 0) then
    do k=1,l
      g = g+z(i,k)*z(i,k)
    end do
  end if
  h = g+f*f
  if (g <= 1d-12) then
    e(i) = f
    h = 0.0d0
  else
    l = l+1
    g = sqrt(h)
    if (f >= 0.0d0) g = -g
    e(i) = g
    h = h-f*g
    z(i,i-1) = f-g
    f = 0.0d0
    do j=1,l
      z(j,i) = z(i,j)/h
      g = 0.0d0
      do k=1,j
        g = g+z(j,k)*z(i,k)
      end do
      jn = j+1
      if (l >= jn) then
        do k=jn,l
          g = g+z(k,j)*z(i,k)
        end do
      end if
      e(j) = g/h
      f = f+g*z(j,i)
    end do
    hh = f/(h+h)
    do j=1,l
      f = z(i,j)
      g = e(j)-hh*f
      e(j) = g
      do k=1,j
        z(j,k) = z(j,k)-f*e(k)-g*z(i,k)
      end do
    end do
  end if
  d(i) = h
end do
d(1) = z(1,1)
z(1,1) = 1.0d0
e(1) = 0.0d0
do i=2,n
  l = i-1
  if (d(i) /= 0.0d0) then
    do j=1,l
      g = 0.0d0
      do k=1,l
        g = g+z(i,k)*z(k,j)
      end do
      do k=1,l
        z(k,j) = z(k,j)-g*z(k,i)
      end do
    end do
  end if
  d(i) = z(i,i)
  z(i,i) = 1.0d0
  do j=1,l
    z(i,j) = 0.0d0
    z(j,i) = 0.0d0
  end do
end do

return

end subroutine hotred

subroutine qlcm(nx,n,d,e,z)

implicit real*8(a-h,o-z)
dimension z(nx,nx), d(nx), e(nx)

do i=2,n
  e(i-1) = e(i)
end do
e(n) = 0.0d0
b = 0.0d0
f = 0.0d0
do l=1,n
  j = 0
  h = 1d-12*(abs(d(l))+abs(e(l)))
  if (b < h) b = h
  do m=l,n
    if (abs(e(m)) <= b) exit
  end do
  if (m /= l) then
    do
      if (j == nx+1) then
        write(6,250)
#       ifndef MOLPRO
        call abend()
#       endif
        !call abend()
      end if
      j = j+1
      g = d(l)
      p = (d(l+1)-g)/(2*e(l))
      r = sqrt(p*p+1.0d0)
      if (p < 0.0d0) then
        pp = p-r
      else
        pp = p+r
      end if
      d(l) = e(l)/pp
      h = g-d(l)
      if (l /= n) then
        ll = l+1
        do i=ll,n
          d(i) = d(i)-h
        end do
      end if
      f = f+h
      p = d(m)
      c = 1.0d0
      s = 0.0d0
      mm = m-1
      do kk=l,mm
        i = mm+l-kk
        g = c*e(i)
        h = c*p
        if (abs(p) < abs(e(i))) then
          c = p/e(i)
          r = sqrt(c*c+1.0d0)
          e(i+1) = s*e(i)*r
          s = 1.0d0/r
          c = c/r
        else
          c = e(i)/p
          r = sqrt(c*c+1.0d0)
          e(i+1) = s*p*r
          s = c/r
          c = 1.0d0/r
        end if
        p = c*d(i)-s*g
        d(i+1) = h+s*(c*g+s*d(i))
        do k=1,n
          h = z(k,i+1)
          z(k,i+1) = s*z(k,i)+c*h
          z(k,i) = c*z(k,i)-s*h
        end do
      end do
      e(l) = s*p
      d(l) = c*p
      if (abs(e(l)) <= b) exit
    end do
  end if
  d(l) = d(l)+f
end do
do i=1,n
  k = i
  p = d(i)
  l = i+1
  if (l <= n) then
    do j=l,n
      if (d(j) >= p) cycle
      k = j
      p = d(j)
    end do
  end if
  if (k == i) cycle
  d(k) = d(i)
  d(i) = p
  do j=1,n
    p = z(j,i)
    z(j,i) = z(j,k)
    z(j,k) = p
  end do
end do

return

250 format(1x///5x,'***the subroutine qlcm is fail, so this computation must stop***')

end subroutine qlcm
