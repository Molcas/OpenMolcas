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

function TestMinimaxLaplace(Tolerance,Verbose)
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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: TestMinimaxLaplace
real(kind=wp), intent(in) :: Tolerance
logical(kind=iwp), intent(in) :: Verbose
integer(kind=iwp) :: irc, K_Lap, l_t, l_t_ref, l_w, l_w_ref, l_wt
real(kind=wp) :: Emax, Emin, RMSt, RMSw, Tol
real(kind=wp), allocatable :: tmlwr(:), tmltr(:), tmlw(:), tmlt(:)
integer(kind=iwp), parameter :: N = 8
real(kind=wp), parameter :: xmax = 1.08976414_wp, xmin = 1.08976414_wp
real(kind=wp), external :: ddot_

if (Verbose) then
  write(u6,'(//,A)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<'
  write(u6,'(A)') '>>>>>>>>>> Enter TestMinimaxLaplace <<<<<<<<<<'
  write(u6,'(A,//)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<'
  call xFlush(u6)
end if

if (Tolerance < Zero) then
  Tol = 1.0e-7_wp
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

tmlwr(1) = 0.097293042801116517_wp
tmlwr(2) = 0.237233954792368945_wp
tmlwr(3) = 0.407050561470959249_wp
tmlwr(4) = 0.635894859784969846_wp
tmlwr(5) = 0.973101625527261427_wp
tmlwr(6) = 1.505487565587009913_wp
tmlwr(7) = 2.419319385961121949_wp
tmlwr(8) = 4.393171503147499379_wp

tmltr(1) = 0.03771106629456916_wp
tmltr(2) = 0.20333950465288628_wp
tmltr(3) = 0.52200685559948912_wp
tmltr(4) = 1.03690003666258002_wp
tmltr(5) = 1.82953857756319826_wp
tmltr(6) = 3.04727454214513793_wp
tmltr(7) = 4.96421466224179220_wp
tmltr(8) = 8.21146012481546705_wp

K_Lap = N
Emin = xmin
Emax = xmax
l_wt = N
call MinimaxLaplace(Verbose,K_Lap,Emin,Emax,l_wt,tmlw,tmlt,irc)
if (Verbose) then
  write(u6,'(A,I6)') 'Return code from MinimaxLaplace=',irc
  call xFlush(u6)
end if
if (irc /= 0) then
  irc = -1
else
  tmlw(1:N) = tmlw(1:N)-tmlwr(1:N)
  tmlt(1:N) = tmlt(1:N)-tmltr(1:N)
  RMSw = sqrt(dDot_(N,tmlw,1,tmlw,1)/real(N,kind=wp))
  RMSt = sqrt(dDot_(N,tmlt,1,tmlt,1)/real(N,kind=wp))
  if (Verbose) then
    write(u6,'(A,1P,D25.16)') 'Weight RMS error=    ',RMSw
    write(u6,'(A,1P,D25.16)') 'Grid point RMS error=',RMSt
    write(u6,'(A,1P,D25.16)') 'Tolerance=           ',Tol
    call xFlush(u6)
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
  write(u6,'(A,I3)') 'TestMinimaxLaplace=',irc
  write(u6,'(//,A)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<'
  write(u6,'(A)') '>>>>>>>>>> Exit TestMinimaxLaplace <<<<<<<<<<'
  write(u6,'(A,//)') '>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<'
  call xFlush(u6)
end if

TestMinimaxLaplace = irc

end function TestMinimaxLaplace
