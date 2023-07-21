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

integer function TestMinimaxLaplace(Tolerance,Verbose)
!
! Thomas Bondo Pedersen, December 2012.
!
! Test that the minimax code from Ten-no gives correct results
! (test case provided by Ten-no). The minimax code is called
! through routine MinimaxLaplace.
!
! Input:
!    Tolerance -- >=0: defines tolerance for test
!                 <0:  default tolerance used
!    Verbose   -- If .true., print information
!
! Returns:
!    -1: nonzero return code from MinimaxLaplace
!     0: weights OK, grid points OK
!     1: weights wrong, grid points OK
!     2: weights OK, grid points wrong
!     3: weights wrong, grid points wrong

use stdalloc

implicit none
real*8 Tolerance
logical Verbose
real*8, allocatable :: tmlwr(:), tmltr(:), tmlw(:), tmlt(:)
real*8 dDot_
external ddot_
integer N
parameter(N=8)
real*8 xmin, xmax
parameter(xmin=1.08976414d0,xmax=1.08976414d0)
integer K_Lap
integer irc
integer l_w_ref
integer l_t_ref
integer l_w
integer l_t
integer l_wt
real*8 Tol
real*8 Emin, Emax
real*8 RMSw, RMSt

if (Verbose) then
  write(6,'(//,A)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<'
  write(6,'(A)') '>>>>>>>>>> Enter TestMinimaxLaplace <<<<<<<<<<'
  write(6,'(A,//)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<'
  call xFlush(6)
end if

if (Tolerance < 0.0d0) then
  Tol = 1.0d-7
else
  Tol = Tolerance
end if

l_w_ref = N
l_t_ref = N
l_w = N
l_t = N
call mma_allocate(tmlwr,l_w_ref,Label='tmlwr')
call mma_allocate(tmltr,l_t_ref,Label='tmltr')
call mma_allocate(tmlw,l_w,Label='tmlw')
call mma_allocate(tmlt,l_t,Label='tmlt')

tmlwr(1) = 0.097293042801116517
tmlwr(2) = 0.237233954792368945
tmlwr(3) = 0.407050561470959249
tmlwr(4) = 0.635894859784969846
tmlwr(5) = 0.973101625527261427
tmlwr(6) = 1.505487565587009913
tmlwr(7) = 2.419319385961121949
tmlwr(8) = 4.393171503147499379

tmltr(1) = 0.03771106629456916
tmltr(2) = 0.20333950465288628
tmltr(3) = 0.52200685559948912
tmltr(4) = 1.03690003666258002
tmltr(5) = 1.82953857756319826
tmltr(6) = 3.04727454214513793
tmltr(7) = 4.96421466224179220
tmltr(8) = 8.21146012481546705

K_Lap = N
Emin = xmin
Emax = xmax
l_wt = N
call MinimaxLaplace(Verbose,K_Lap,Emin,Emax,l_wt,tmlw,tmlt,irc)
if (Verbose) then
  write(6,'(A,I6)') 'Return code from MinimaxLaplace=',irc
  call xFlush(6)
end if
if (irc /= 0) then
  irc = -1
else
  call dAXPY_(N,-1.0d0,tmlwr,1,tmlw,1)
  call dAXPY_(N,-1.0d0,tmltr,1,tmlt,1)
  RMSw = sqrt(dDot_(N,tmlw,1,tmlw,1)/dble(N))
  RMSt = sqrt(dDot_(N,tmlt,1,tmlt,1)/dble(N))
  if (Verbose) then
    write(6,'(A,1P,D25.16)') 'Weight RMS error=    ',RMSw
    write(6,'(A,1P,D25.16)') 'Grid point RMS error=',RMSt
    write(6,'(A,1P,D25.16)') 'Tolerance=           ',Tol
    call xFlush(6)
  end if
  irc = 0
  if (RMSw > Tol) irc = irc+1
  if (RMSt > Tol) irc = irc+2
end if

call mma_deallocate(tmlt)
call mma_deallocate(tmlw)
call mma_deallocate(tmltr)
call mma_deallocate(tmlwr)

if (Verbose) then
  write(6,'(A,I3)') 'TestMinimaxLaplace=',irc
  write(6,'(//,A)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<'
  write(6,'(A)') '>>>>>>>>>> Exit TestMinimaxLaplace <<<<<<<<<<'
  write(6,'(A,//)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<'
  call xFlush(6)
end if

TestMinimaxLaplace = irc

end function TestMinimaxLaplace
