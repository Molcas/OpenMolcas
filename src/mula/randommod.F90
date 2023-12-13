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

module RandomMod

! Description:
!   This random number generator originally appeared in "Toward a Universal
!   Random Number Generator" by George Marsaglia and Arif Zaman.
!   Florida State University Report: FSU-SCRI-87-50 (1987)
!
!   It was later modified by F. James and published in "A Review of Pseudo-
!   random Number Generators"
!
!   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
!         (However, a newly discovered technique can yield
!           a period of 10^600. But that is still in the development stage.)
!
!   It passes ALL of the tests for random number generators and has a period
!   of 2^144, is completely portable (gives bit identical results on all
!   machines with at least 24-bit mantissas in the floating point
!   representation).
!
!   The algorithm is a combination of a Fibonacci sequence (with lags of 97
!   and 33, and operation "subtraction plus one, modulo one") and an
!   "arithmetic sequence" (using subtraction).
!
!   On a Vax 11/780, this random number generator can produce a number in
!   13 microseconds.

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp) :: i97, j97
real(kind=wp) :: c, cd, cm, U(97)
logical(kind=iwp) :: test = .false.

public :: Ranmar, Rmarin

contains

subroutine Rmarin(ij,kl)
! This is the initialization routine for the random number generator RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
! The random number sequences created by these two seeds are of sufficient
! length to complete an entire calculation with. For example, if sveral
! different groups are working on different parts of the same calculation,
! each group could be assigned its own IJ seed. This would leave each group
! with 30000 choices for the second seed. That is to say, this random
! number generator can create 900 million different subsequences -- with
! each subsequence having a length of approximately 10^30.
!
! Use IJ = 1802 & KL = 9373 to test the random number generator. The
! subroutine RANMAR should be used to generate 20000 random numbers.
! Then display the next six random numbers generated multiplied by 4096*4096
! If the random number generator is working properly, the random numbers
! should be:
!           6533892.0  14220222.0  7275067.0
!           6172232.0  8354498.0   10633180.0

  use Constants, only: Zero, Half
  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: ij, kl
  integer(kind=iwp) :: i, ii, j, jj, k, l, m
  real(kind=wp) :: s, t

  test = .false.

  if ((ij < 0) .or. (ij > 31328) .or. (kl < 0) .or. (kl > 30081)) then
    write(u6,'(a)') 'The first random number seed must have a value between 0 and 31328'
    write(u6,'(a)') 'The second seed must have a value between 0 and 30081'
    call abend()
  end if

  i = mod(ij/177,177)+2
  j = mod(ij,177)+2
  k = mod(kl/169,178)+1
  l = mod(kl,169)

  do ii=1,97
    s = Zero
    t = Half
    do jj=1,24
      m = mod(mod(i*j,179)*k,179)
      i = j
      j = k
      k = m
      l = mod(53*l+1,169)
      if (mod(l*m,64) >= 32) then
        s = s+t
      end if
      t = Half*t
    end do
    U(ii) = s
  end do

  c = 362436.0_wp/16777216.0_wp
  cd = 7654321.0_wp/16777216.0_wp
  cm = 16777213.0_wp/16777216.0_wp

  i97 = 97
  j97 = 33

  test = .true.

end subroutine Rmarin

subroutine Ranmar(rvec,len)
! This is the random number generator proposed by George Marsaglia in
! Florida State University Report: FSU-SCRI-87-50
! It was slightly modified by F. James to produce an array of pseudorandom
! numbers.

  use Constants, only: Zero, One
  use Definitions, only: u6

  integer(kind=iwp), intent(in) :: len
  real(kind=wp), intent(out) :: rvec(len)
  integer(kind=iwp) :: ivec
  real(kind=wp) :: uni

  if (.not. test) then
    write(u6,'(a)') ' Call the init routine (RMARIN) before calling RANMAR'
    call abend()
  end if

  do ivec=1,len
    uni = u(i97)-u(j97)
    if (uni < Zero) uni = uni+One
    u(i97) = uni
    i97 = i97-1
    if (i97 == 0) i97 = 97
    j97 = j97-1
    if (j97 == 0) j97 = 97
    c = c-cd
    if (c < Zero) c = c+cm
    uni = uni-c
    if (uni < Zero) uni = uni+One
    rvec(ivec) = uni
  end do

end subroutine Ranmar

end module RandomMod
