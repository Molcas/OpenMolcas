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
!     labs(maxlab)           labels for components of the l-pole
!                            moment
!     ndim                   defines the number of rows which will
!                            be transformed in the table
!                            prvec(ndim,maxlab)
!     t(maxlab,maxlab)       is the transformation matrix generated
!                            in this program for lpole=2,3,4
!     temp(maxlab)           is a temporary strorage area
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, Two, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: inp, lpole, maxlab, ndim
character(len=16), intent(in) :: labs(maxlab)
real(kind=wp), intent(inout) :: prvec(ndim,maxlab), t(maxlab,maxlab)
real(kind=wp), intent(out) :: temp(maxlab)
integer(kind=iwp) :: i, icount, ilab(6,3), ind, irr(3,3), irrrr(6,3), j, k, l
real(kind=wp) :: f, rsum
character :: l1, l2, l3, l4

k = 0 ! dummy initialize
l = 0 ! dummy initialize

! building of transformation matrices for the transformation
! from cartesian l-th moments to cartesian l-pole moments
! limited to l=2,3, and 4

if (inp /= 1) then

  select case (lpole)

    case (2)
      ! quadrupole moments

      t(:,:) = Zero
      do i=1,maxlab
        t(i,i) = OneHalf
        read(labs(i),'(14x,2a1)') l1,l2
        if (l1 == l2) then
          t(i,1) = t(i,1)-Half
          t(i,4) = t(i,4)-Half
          t(i,6) = t(i,6)-Half
        end if
      end do

    case (3)
      ! octupole moments

      irr(:,:) = 0
      do i=1,3
        irr(i,i) = 2
      end do

      t(:,:) = Zero
      do i=1,maxlab
        t(i,i) = 2.5_wp
        read(labs(i),'(13x,3a1)') l1,l2,l3
        if (l1 == l2) then
          ilab(1:3,:) = irr
          k = 0
          if (l3 == 'X') k = 1
          if (l3 == 'Y') k = 2
          if (l3 == 'Z') k = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          do j=1,3
            ind = nTri_Elem(3-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-Half
          end do
        end if
        if (l2 == l3) then
          ilab(1:3,:) = irr
          k = 0
          if (l1 == 'X') k = 1
          if (l1 == 'Y') k = 2
          if (l1 == 'Z') k = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          do j=1,3
            ind = nTri_Elem(3-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-Half
          end do
        end if
        if (l1 == l3) then
          ilab(1:3,:) = irr
          k = 0
          if (l2 == 'X') k = 1
          if (l2 == 'Y') k = 2
          if (l2 == 'Z') k = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          do j=1,3
            ind = nTri_Elem(3-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-Half
          end do
        end if
      end do

    case (4)
      ! hexadecapole moments

      irr(:,:) = 0
      do i=1,3
        irr(i,i) = 2
      end do
      irrrr(:,:) = 0
      irrrr(1,1) = 4
      irrrr(2,1) = 2
      irrrr(2,2) = 2
      irrrr(3,1) = 2
      irrrr(3,3) = 2
      irrrr(4,2) = 4
      irrrr(5,2) = 2
      irrrr(5,3) = 2
      irrrr(6,3) = 4

      t(:,:) = Zero
      do i=1,maxlab
        t(i,i) = 4.375_wp
        read(labs(i),'(12x,4a1)') l1,l2,l3,l4
        if (l1 == l2) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l3 == 'X') k = 1
          if (l3 == 'Y') k = 2
          if (l3 == 'Z') k = 3
          if (l4 == 'X') l = 1
          if (l4 == 'Y') l = 2
          if (l4 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        if (l1 == l3) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l2 == 'X') k = 1
          if (l2 == 'Y') k = 2
          if (l2 == 'Z') k = 3
          if (l4 == 'X') l = 1
          if (l4 == 'Y') l = 2
          if (l4 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        if (l1 == l4) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l2 == 'X') k = 1
          if (l2 == 'Y') k = 2
          if (l2 == 'Z') k = 3
          if (l3 == 'X') l = 1
          if (l3 == 'Y') l = 2
          if (l3 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        if (l2 == l3) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l1 == 'X') k = 1
          if (l1 == 'Y') k = 2
          if (l1 == 'Z') k = 3
          if (l4 == 'X') l = 1
          if (l4 == 'Y') l = 2
          if (l4 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        if (l2 == l4) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l1 == 'X') k = 1
          if (l1 == 'Y') k = 2
          if (l1 == 'Z') k = 3
          if (l3 == 'X') l = 1
          if (l3 == 'Y') l = 2
          if (l3 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        if (l3 == l4) then
          ilab(1:3,:) = irr
          k = 0
          l = 0
          if (l1 == 'X') k = 1
          if (l1 == 'Y') k = 2
          if (l1 == 'Z') k = 3
          if (l2 == 'X') l = 1
          if (l2 == 'Y') l = 2
          if (l2 == 'Z') l = 3
          ilab(1:3,k) = ilab(1:3,k)+1
          ilab(1:3,l) = ilab(1:3,l)+1
          do j=1,3
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)-0.625_wp
          end do
        end if
        ilab(:,:) = irrrr
        if ((l1 == l2) .and. (l3 == l4)) then
          do j=1,6
            f = 0.125_wp
            if ((j == 2) .or. (j == 3) .or. (j == 5)) f = Two*f
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)+f
          end do
        end if
        if ((l1 == l3) .and. (l2 == l4)) then
          do j=1,6
            f = 0.125_wp
            if ((j == 2) .or. (j == 3) .or. (j == 5)) f = Two*f
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)+f
          end do
        end if
        if ((l2 == l3) .and. (l1 == l4)) then
          do j=1,6
            f = 0.125_wp
            if ((j == 2) .or. (j == 3) .or. (j == 5)) f = Two*f
            ind = nTri_Elem(4-ilab(j,1))+ilab(j,3)+1
            T(i,ind) = T(i,ind)+f
          end do
        end if
      end do

    case default
      call Abend()

  end select

  ! print the transformation matrix

  !write(u6,'(//1x,a,i2/)') 'transformation matrix:  lpole=',lpole
  !do i=1,maxlab
  !  write(u6,'(15f7.3)') t(i,:)
  !end do

end if

! transform cartesian moment to multipole moments

do icount=1,ndim
  temp(:) = prvec(icount,:)
  do k=1,maxlab
    rsum = Zero
    do l=1,maxlab
      rsum = rsum+t(k,l)*temp(l)
    end do
    prvec(icount,k) = rsum
  end do
end do

return

end subroutine tmltpl
