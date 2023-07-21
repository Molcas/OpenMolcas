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
! Copyright (C) 2012, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine MinimaxLaplace(Verbose,N,xmin,xmax,l_wt,w,t,irc)
!
! Thomas Bondo Pedersen, November 2012.
!
! Generate grid (weights w and points t) for the
! Laplace transform of 1/x.
!
! The minimax approximation for setting up the grid
! [Takatsuka, Ten-no, Hackbusch; JCP 129, 044112 (2008)]
! is used. The source code for the minimax approximation
! was provided by Ten-no on Nov. 28, 2012 (to TBP).
!
! Verbose   -- boolean to control printing
! N         -- number of grid points requested. If N=0 on input,
!              the number of points required to get an accuracy
!              of 1.0d-6 (in the Laplace transform) is used. In
!              this case, N is the required number of grid points
!              on return.
! xmin,xmax -- range of x values to be covered.
! l_wt      -- allocated dimension of the arrays w and t.
!              If N(input)=0 and l_wt<N(output), only l_wt grid
!              points are returned in w and t. N(output) is still
!              the full number of grid points [this situation is
!              also flagged by the return code irc, see below].
! w         -- grid weights are returned in array w.
! t         -- grid points are returned in array t.
! irc       -- return code:
!              <0: input error (parameter -irc is the culprit)
!               0: all OK
!               1: l_wt is too small to store full result and only
!                  l_wt grid weights/points are returned.

use stdalloc

implicit none
logical Verbose
integer N
real*8 xmin
real*8 xmax
integer l_wt
real*8 w(l_wt)
real*8 t(l_wt)
integer irc
integer mGrid
parameter(mGrid=20) ! limited by Remez implementation
character*8 DefaultGrid
parameter(DefaultGrid='MICRO   ')
integer l_Coeff
real*8, allocatable :: Coeff(:)
integer K_Lap
character*8 Demand
logical Inf
integer i

irc = 0
if ((N < 0) .or. (N > mGrid)) then
  irc = -1
  return
end if
if (xmin < 0.0d0) then
  irc = -2
  return
end if
if ((xmax-xmin) < 0.0d0) then
  irc = -3
  return
end if
if (l_wt < 1) then
  irc = -4
  return
end if

K_Lap = N
if (N == 0) then
  Demand = DefaultGrid
else
  Demand = '        '
end if
l_Coeff = 2*mGrid
call mma_Allocate(Coeff,l_Coeff,Label='LapCoef')
Inf = .false.
call Remez(Verbose,K_Lap,xmin,xmax,Coeff,Demand,Inf)
if (K_Lap < 0) then
  call mma_Deallocate(Coeff)
  irc = -1
  write(6,'(A,I10)') 'MinimaxLaplace: Remez returned K_Lap=',K_Lap
  return
end if
if (N == 0) N = K_Lap
if (l_wt < K_Lap) then
  do i=1,l_wt
    w(i) = Coeff(2*i-1)
    t(i) = Coeff(2*i)
  end do
  irc = 2
else
  do i=1,K_Lap
    w(i) = Coeff(2*i-1)
    t(i) = Coeff(2*i)
  end do
end if
call mma_Deallocate(Coeff)

end subroutine MinimaxLaplace
