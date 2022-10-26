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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine tmltpl(inp,lpole,maxlab,labs,ndim,prvec,t,temp)
!***********************************************************************
!
!     Purpose: transformation of maxlab cartesian l-th
!              moments into the corresponding l-pole cartesian
!              moments; l=lpole
!
!
!     inp                    if inp == 0 the t matrix will be
!                            calculated. Otherwise the t matrix
!                            is transferred from the previous
!                            call
!     lpole                  l value for the l-pole moment
!     maxlab                 is the number of cartesian components
!                            of the l-pole moment
!     labs(1:maxlab)         labels for components of the l-pole
!                            moment
!     ndim                   defines the number of rows which will
!                            be transformed in the table
!                            prvec(1:ndim,1:maxlab)
!     t(1:maxlab,1:maxlab)   is the transformation matrix generated
!                            in this program for lpole=2,3,4
!     temp(1:maxlab)         is a temporary strorage area
!
!***********************************************************************

implicit real*8(a-h,o-z)
character*1 l1, l2, l3, l4
character*14 l14
character*16 labs(1:maxlab)
dimension prvec(1:ndim,1:maxlab)
dimension t(1:maxlab,1:maxlab)
dimension temp(1:maxlab)
dimension irr(1:3,1:3), ilab(1:6,1:3), irrrr(1:6,1:3)

k = 0 ! dummy initialize
l = 0 ! dummy initialize

! building of transformation matrices for the transformation
! from cartesian l-th moments to cartesian l-pole moments
! limited to l=2,3, and 4

if (inp == 1) go to 98
go to(100,200,300),lpole-1

! quadrupole moments

100 continue

do i=1,maxlab
  do j=1,maxlab
    t(i,j) = 0.0d+00
  end do
  t(i,i) = t(i,i)+1.5d+00
  read(labs(i),'(a14,2a1)') l14,l1,l2
  if (l1 == l2) then
    t(i,1) = t(i,1)-0.5d+00
    t(i,4) = t(i,4)-0.5d+00
    t(i,6) = t(i,6)-0.5d+00
  end if
end do
go to 99

! octupole moments

200 continue

do i=1,3
  do j=1,3
    irr(i,j) = 0
  end do
  irr(i,i) = 2
end do

do i=1,maxlab
  do j=1,maxlab
    t(i,j) = 0.0d+00
  end do
  t(i,i) = t(i,i)+2.5d+00
  read(labs(i),'(a13,3a1)') l14,l1,l2,l3
  if (l1 == l2) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    if (l3 == 'X') k = 1
    if (l3 == 'Y') k = 2
    if (l3 == 'Z') k = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ind = (3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.5d+00
    end do
  end if
  if (l2 == l3) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    if (l1 == 'X') k = 1
    if (l1 == 'Y') k = 2
    if (l1 == 'Z') k = 3
    do j=1,3
      ilab(j,k) = irr(j,k)+1
      ind = (3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.5d+00
    end do
  end if
  if (l1 == l3) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    if (l2 == 'X') k = 1
    if (l2 == 'Y') k = 2
    if (l2 == 'Z') k = 3
    do j=1,3
      ilab(j,k) = irr(j,k)+1
      ind = (3-ilab(j,1))*(3-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.5d+00
    end do
  end if
end do
go to 99

! hexadecapole moments

300 continue

do i=1,3
  do j=1,3
    irr(i,j) = 0
  end do
  do j=1,6
    irrrr(j,i) = 0
  end do
  irr(i,i) = 2
end do
irrrr(1,1) = 4
irrrr(2,1) = 2
irrrr(2,2) = 2
irrrr(3,1) = 2
irrrr(3,3) = 2
irrrr(4,2) = 4
irrrr(5,2) = 2
irrrr(5,3) = 2
irrrr(6,3) = 4

do i=1,maxlab
  do j=1,maxlab
    t(i,j) = 0.0d+00
  end do
  t(i,i) = t(i,i)+4.375d+00
  read(labs(i),'(a12,4a1)') l14,l1,l2,l3,l4
  if (l1 == l2) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l3 == 'X') k = 1
    if (l3 == 'Y') k = 2
    if (l3 == 'Z') k = 3
    if (l4 == 'X') l = 1
    if (l4 == 'Y') l = 2
    if (l4 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  if (l1 == l3) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l2 == 'X') k = 1
    if (l2 == 'Y') k = 2
    if (l2 == 'Z') k = 3
    if (l4 == 'X') l = 1
    if (l4 == 'Y') l = 2
    if (l4 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  if (l1 == l4) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l2 == 'X') k = 1
    if (l2 == 'Y') k = 2
    if (l2 == 'Z') k = 3
    if (l3 == 'X') l = 1
    if (l3 == 'Y') l = 2
    if (l3 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  if (l2 == l3) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l1 == 'X') k = 1
    if (l1 == 'Y') k = 2
    if (l1 == 'Z') k = 3
    if (l4 == 'X') l = 1
    if (l4 == 'Y') l = 2
    if (l4 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  if (l2 == l4) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l1 == 'X') k = 1
    if (l1 == 'Y') k = 2
    if (l1 == 'Z') k = 3
    if (l3 == 'X') l = 1
    if (l3 == 'Y') l = 2
    if (l3 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  if (l3 == l4) then
    do i1=1,3
      do i2=1,3
        ilab(i1,i2) = irr(i1,i2)
      end do
    end do
    k = 0
    l = 0
    if (l1 == 'X') k = 1
    if (l1 == 'Y') k = 2
    if (l1 == 'Z') k = 3
    if (l2 == 'X') l = 1
    if (l2 == 'Y') l = 2
    if (l2 == 'Z') l = 3
    do j=1,3
      ilab(j,k) = ilab(j,k)+1
      ilab(j,l) = ilab(j,l)+1
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)-0.625d+00
    end do
  end if
  do i1=1,6
    do i2=1,3
      ilab(i1,i2) = irrrr(i1,i2)
    end do
  end do
  if ((l1 == l2) .and. (l3 == l4)) then
    do j=1,6
      f = 0.125d+00
      if ((j == 2) .or. (j == 3) .or. (j == 5)) f = 2.0d+00*f
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)+f
    end do
  end if
  if ((l1 == l3) .and. (l2 == l4)) then
    do j=1,6
      f = 0.125d+00
      if ((j == 2) .or. (j == 3) .or. (j == 5)) f = 2.0d+00*f
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)+f
    end do
  end if
  if ((l2 == l3) .and. (l1 == l4)) then
    do j=1,6
      f = 0.125d+00
      if ((j == 2) .or. (j == 3) .or. (j == 5)) f = 2.0d+00*f
      ind = (4-ilab(j,1))*(4-ilab(j,1)+1)/2+ilab(j,3)+1
      T(i,ind) = T(i,ind)+f
    end do
  end if
end do

99 continue

! print the transformation matrix

!write(6,'(//1x,a,i2/)') 'transformation matrix:  lpole=',lpole
!do i=1,maxlab
!  write(6,'(15f7.3)') (t(i,j),j=1,maxlab)
!end do

98 continue

! transform cartesian moment to multipole moments

do icount=1,ndim
  do k=1,maxlab
    temp(k) = prvec(icount,k)
  end do
  do k=1,maxlab
    sum = +0.0d+00
    do l=1,maxlab
      sum = sum+t(k,l)*temp(l)
    end do
    prvec(icount,k) = sum
  end do
end do

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call unused_character(l14)
#endif

end subroutine tmltpl
