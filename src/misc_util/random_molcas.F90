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
! Copyright (C) 1993, Per-Olof Widmark                                 *
!               2010, John Burkardt                                    *
!               2017, Morgane Vacher                                   *
!***********************************************************************

function Random_Molcas(iSeed)
!***********************************************************************
!                                                                      *
!     Generate a random number z, where 0<z<0.5**31                    *
!                                                                      *
!     calling arguments:                                               *
!     iSeed  : Type integer, input/output                              *
!              initial value of new chain                              *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     P.O. Widmark                                                     *
!     University of Lund, Sweden, 1993                                 *
!     M. Vacher, Uppsala 2017                                          *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use Definitions, only: wp, iwp
use, intrinsic :: iso_fortran_env, only: real64

implicit none
real(kind=wp) :: Random_Molcas
integer(kind=iwp), intent(inout) :: iSeed
integer(kind=iwp) :: i, IX0, IX1, IX2, IX3
real(kind=real64) :: a1 = 0.0_real64, a2 = 0.0_real64, r23 = 1.0_real64, r46 = 1.0_real64, t1, t2, t23 = 1.0_real64, t3, t4, &
                     t46 = 1.0_real64, x, x1, x2, z
logical(kind=iwp) :: ks = .true.
character(len=8) :: rand_info
integer(kind=iwp), parameter :: IA1 = 8121, IA2 = 4561, IA3 = 7141, IC1 = 28411, IC2 = 51349, IC3 = 54773, IM1 = 134456, &
                                IM2 = 243000, IM3 = 259200
real(kind=real64), parameter :: a = 1220703125.0_real64, x0 = 314159265.0_real64

call GetEnvf('MOLCAS_RANDOM',rand_info)
call UpCase(rand_info)
if (rand_info(1:3) == 'OLD') then
  IX0 = iSeed
  IX1 = mod(IA1*IX0+IC1,IM1)
  IX2 = mod(IA2*IX1+IC2,IM2)
  IX3 = mod(IA3*IX2+IC3,IM3)
  Random_Molcas = (real(IX1,kind=wp)+real(IX2,kind=wp)/real(IM2,kind=wp))/real(IM1,kind=wp)
  iSeed = IX3
else
  !*********************************************************************
  ! RANDLC                                                             *
  ! The number returned is a uniform pseudorandom value in the range   *
  ! (0,1). The algorithm uses the linear congruential generator:       *
  !   X(K+1) = A * X(K)  mod 2^46                                      *
  ! This scheme generates 2^44 numbers before repeating.               *
  !--------------------------------------------------------------------*
  ! Author: John Burkardt                                              *
  ! Modified: 08 March 2010                                            *
  !--------------------------------------------------------------------*
  ! References:                                                        *
  ! * David Bailey, Eric Barszcz, John Barton, D Browning, Robert      *
  ! Carter, Leonardo Dagum, Rod Fatoohi, Samuel Fineberg, Paul         *
  ! Frederickson, Thomas Lasinski, Robert Schreiber, Horst Simon, V    *
  ! Venkatakrishnan, Sisira Weeratunga, The NAS Parallel Benchmarks,   *
  ! RNR Technical Report RNR-94-007, March 1994.                       *
  ! * Donald Knuth, The Art of Computer Programming, Volume 2,         *
  ! Seminumerical Algorithms, Third Edition, Addison Wesley, 1997,     *
  ! ISBN: 0201896842, LC: QA76.6.K64.                                  *
  !--------------------------------------------------------------------*
  ! Parameters:                                                        *
  ! Input/output, real(kind=real64) X, the seed.  X should be an odd   *
  ! integer such that 1 <= X <= 2^46.                                  *
  ! Output, real(kind=wp), the next pseudorandom value.                *
  !*********************************************************************
  ! Should we ensure X is odd?
  x = real(iseed,kind=real64)
  ! If this is the first call, compute
  !   R23 = 2 ^ -23,
  !   R46 = 2 ^ -46,
  !   T23 = 2 ^ 23,
  !   T46 = 2 ^ 46.
  ! These are computed in loops, rather than by merely using the power operator,
  ! in order to insure that the results are exact on all systems.
  if (ks) then
    do i=1,46
      r46 = 0.5_real64*r46
      t46 = 2.0_real64*t46
      if (i == 23) then
        r23 = r46
        t23 = t46
      end if
    end do
    ! Break A into two parts such that A = 2^23 * A1 + A2.
    t1 = r23*a
    a1 = real(int(t1),kind=real64)
    a2 = a-t23*a1
    ks = .false.
  end if
  ! Deal with a 0 input value of X.
  if (x == 0.0_real64) x = x0
  ! Deal somewhat arbitrarily with negative input X.
  if (x < 0.0_real64) x = -x
  ! Break X into two parts X1 and X2 such that:
  !   X = 2^23 * X1 + X2,
  ! then compute
  !   Z = A1 * X2 + A2 * X1  (mod 2^23)
  !   X = 2^23 * Z + A2 * X2  (mod 2^46).
  t1 = r23*x
  x1 = real(int(t1),kind=real64)
  x2 = x-t23*x1
  t1 = a1*x2+a2*x1
  t2 = real(int(r23*t1),kind=real64)
  z = t1-t23*t2
  t3 = t23*z+a2*x2
  t4 = real(int(r46*t3),kind=real64)
  x = t3-t46*t4
  Random_Molcas = real(r46*x,kind=wp)
  ! Should we ensure iseed is odd?
  iseed = int(x,kind=iwp)
end if

!----------------------------------------------------------------------*
! exit                                                                 *
!----------------------------------------------------------------------*
return

end function Random_Molcas
