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

module wigner_util
! Thin OpenMolcas-side wrappers around libwignernj, which evaluates the
! Wigner 3-j/6-j/9-j symbols and Clebsch-Gordan coefficients exactly
! (prime-factorised, multi-word integer Racah summation).
!
! libwignernj takes every angular-momentum argument as an integer equal
! to twice its true (possibly half-integer) value -- the "2*j" integer
! convention -- through a C-interoperable interface. These wrappers
! convert from OpenMolcas's call-site conventions so that callers deal
! with neither c_int kinds nor the doubling. Each symbol exposes
! dedicated names for the real-valued and 2*j integer forms rather
! than a single generic interface, since at least one compiler (NVHPC)
! mis-dispatches such generic calls under -i8.
!
!   w3j     -- Wigner 3-j symbol ( j1 j2 j3 / m1 m2 m3 ), real arguments
!              are the actual (half-)integer momenta
!   w3j_2j  -- same, integer arguments already in the 2*j convention
!   w6j     -- Wigner 6-j symbol { j1 j2 j3 / j4 j5 j6 }, integer 2*j
!   w9j     -- Wigner 9-j symbol ( 3x3 array of momenta ), integer 2*j
!   wcg     -- Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>,
!              arguments interleaved (each m next to its j), integer 2*j
!   wcg_real-- same, real arguments are actual (half-)integer momenta
!   dclebs  -- Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>,
!              arguments grouped as (j1,j2,j3,m1,m2,m3), real

use, intrinsic :: iso_c_binding, only: c_int
use wignernj, only: clebsch_gordan, wigner3j, wigner6j, wigner9j
use Constants, only: Two
use Definitions, only: wp, iwp

implicit none
private

public :: dclebs, w3j, w3j_2j, w6j, w9j, wcg, wcg_real

contains

function w3j(j1,j2,j3,m1,m2,m3)
  real(kind=wp) :: w3j
  real(kind=wp), intent(in) :: j1, j2, j3, m1, m2, m3
  w3j = wigner3j(nint(Two*j1,c_int),nint(Two*j2,c_int),nint(Two*j3,c_int), &
                 nint(Two*m1,c_int),nint(Two*m2,c_int),nint(Two*m3,c_int))
end function w3j

function w3j_2j(j1,j2,j3,m1,m2,m3)
  ! Arguments are in the 2*j convention.
  real(kind=wp) :: w3j_2j
  integer(kind=iwp), intent(in) :: j1, j2, j3, m1, m2, m3
  w3j_2j = wigner3j(int(j1,c_int),int(j2,c_int),int(j3,c_int), &
                    int(m1,c_int),int(m2,c_int),int(m3,c_int))
end function w3j_2j

function w6j(j1,j2,j3,j4,j5,j6)
  ! Arguments are in the 2*j convention.
  real(kind=wp) :: w6j
  integer(kind=iwp), intent(in) :: j1, j2, j3, j4, j5, j6
  w6j = wigner6j(int(j1,c_int),int(j2,c_int),int(j3,c_int), &
                 int(j4,c_int),int(j5,c_int),int(j6,c_int))
end function w6j

function w9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
  ! Arguments are in the 2*j convention.
  real(kind=wp) :: w9j
  integer(kind=iwp), intent(in) :: j1, j2, j3, j4, j5, j6, j7, j8, j9
  w9j = wigner9j(int(j1,c_int),int(j2,c_int),int(j3,c_int), &
                 int(j4,c_int),int(j5,c_int),int(j6,c_int), &
                 int(j7,c_int),int(j8,c_int),int(j9,c_int))
end function w9j

function wcg(j1,m1,j2,m2,j3,m3)
  ! Arguments are in the 2*j convention.
  real(kind=wp) :: wcg
  integer(kind=iwp), intent(in) :: j1, m1, j2, m2, j3, m3
  wcg = clebsch_gordan(int(j1,c_int),int(m1,c_int),int(j2,c_int), &
                       int(m2,c_int),int(j3,c_int),int(m3,c_int))
end function wcg

function wcg_real(j1,m1,j2,m2,j3,m3)
  real(kind=wp) :: wcg_real
  real(kind=wp), intent(in) :: j1, m1, j2, m2, j3, m3
  wcg_real = clebsch_gordan(nint(Two*j1,c_int),nint(Two*m1,c_int), &
                            nint(Two*j2,c_int),nint(Two*m2,c_int), &
                            nint(Two*j3,c_int),nint(Two*m3,c_int))
end function wcg_real

function dclebs(j1,j2,j3,m1,m2,m3)
  ! Clebsch-Gordan coefficient <j1 m1 j2 m2 | j3 m3>, arguments grouped.
  real(kind=wp) :: dclebs
  real(kind=wp), intent(in) :: j1, j2, j3, m1, m2, m3
  dclebs = clebsch_gordan(nint(Two*j1,c_int),nint(Two*m1,c_int), &
                          nint(Two*j2,c_int),nint(Two*m2,c_int), &
                          nint(Two*j3,c_int),nint(Two*m3,c_int))
end function dclebs

end module wigner_util
