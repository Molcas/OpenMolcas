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

subroutine DGETMO(A,ldA,M,N,B,ldB)
! TRANSPOSE A REGULAR MATRIX (OUT-OF-PLACE)

#include "intent.fh"

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ldA, M, N, ldB
real(kind=wp), intent(in) :: A(ldA,N)
real(kind=wp), intent(_OUT_) :: B(ldB,M)
integer(kind=iwp) :: i, INC, j, jj

if (M <= 0) then
  write(u6,*)
  write(u6,*) '  *** Error in subroutine DGETMO ***'
  write(u6,*) '  Invalid dimension of matrix A :'
  write(u6,*) '  The number of columns, M, must be greater than zero'
  write(u6,*)
  call abend()
end if
if (N <= 0) then
  write(u6,*)
  write(u6,*) '  *** Error in subroutine DGETMO ***'
  write(u6,*) '  Invalid leading dimension of matrix B :'
  write(u6,*) '  The number of rows, N, must be greater than zero'
  write(u6,*)
  call abend()
end if
if (ldA < M) then
  write(u6,*)
  write(u6,*) '  *** Error in subroutine DGETMO ***'
  write(u6,*) '  Invalid leading dimension of matrix A :'
  write(u6,*) '  ldA must be equal to M or greater'
  write(u6,*)
  call abend()
end if
if (ldB < N) then
  write(u6,*)
  write(u6,*) '  *** Error in subroutine DGETMO ***'
  write(u6,*) '  Invalid leading dimension of matrix B :'
  write(u6,*) '  ldB must be equal to N or greater'
  write(u6,*)
  call abend()
end if
INC = 8
do j=1,M,INC
  jj = min(M-j+1,INC)
  select case (jj)
    case (1)
      do i=1,N
        B(i,j) = A(j,i)
      end do
    case (2)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
      end do
    case (3)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
      end do
    case (4)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
        B(i,j+3) = A(j+3,i)
      end do
    case (5)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
        B(i,j+3) = A(j+3,i)
        B(i,j+4) = A(j+4,i)
      end do
    case (6)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
        B(i,j+3) = A(j+3,i)
        B(i,j+4) = A(j+4,i)
        B(i,j+5) = A(j+5,i)
      end do
    case (7)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
        B(i,j+3) = A(j+3,i)
        B(i,j+4) = A(j+4,i)
        B(i,j+5) = A(j+5,i)
        B(i,j+6) = A(j+6,i)
      end do
    case (8)
      do i=1,N
        B(i,j) = A(j,i)
        B(i,j+1) = A(j+1,i)
        B(i,j+2) = A(j+2,i)
        B(i,j+3) = A(j+3,i)
        B(i,j+4) = A(j+4,i)
        B(i,j+5) = A(j+5,i)
        B(i,j+6) = A(j+6,i)
        B(i,j+7) = A(j+7,i)
      end do
    case default
      write(u6,*) 'Error in DGETMO!'
  end select
end do

return

end subroutine DGETMO
