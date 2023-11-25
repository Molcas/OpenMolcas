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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine tunedefs2_cvb(imethod,endwhenclose1)

use casvb_global, only: delopth1, delopth2, dfxmin, dx, endwhenclose, exp12tol, grd, hhaccfac, hhmax, hhrejfac, hhstart, hhtol, &
                        nopth1, nopth2, resthr, scalesmall, sgn, singul, zzacclim, zzmin, zzrejmax, zzrejmin

use Constants, only: Zero, One, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: imethod
logical(kind=iwp), intent(in) :: endwhenclose1
real(kind=wp), parameter :: small = 1.0e-3_wp, small2 = 1.0e-5_wp, smaller = 1.0e-6_wp

! << TUNE_CVB common block: >>
endwhenclose = endwhenclose1
if ((imethod == 1) .or. (imethod == 10)) then
  !'FLETCHER'
  exp12tol = Zero
  ! Criteria for entering local region:
  grd(1,1) = 5.0e-4_wp
  grd(1,2) = 5.0e-4_wp
  sgn(1) = smaller
  sgn(2) = smaller
  singul(1) = 1.0e-2_wp
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  grd(1,3) = 5.0e-6_wp
  grd(1,4) = 5.0e-6_wp
  dx(1,3) = 5.0e-6_wp
  dx(1,4) = 1.0e-4_wp
  singul(2) = 1.0e-2_wp
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = Zero
  hhrejfac(1) = 0.4_wp
  hhaccfac(2,1) = One
  hhaccfac(3,1) = OneHalf
  hhaccfac(4,1) = One
  zzacclim(2,2) = 0.8_wp
  zzacclim(3,2) = 1.25_wp
  hhtol(1) = 1.0e-10_wp
  dfxmin(1) = Zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = 0.4_wp
  hhaccfac(3,2) = 1.2_wp
  hhtol(2) = 1.0e-10_wp
  dfxmin(2) = Zero
end if
if (imethod == 2) then
  ! 'TRIM'
  exp12tol = Zero
  ! Criteria for entering local region:
  grd(1,1) = 5.0e-6_wp
  grd(1,2) = 5.0e-6_wp
  sgn(1) = smaller
  sgn(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  dx(1,3) = 5.0e-6_wp
  dx(1,4) = 1.0e-4_wp
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = 0.75_wp
  zzrejmax(1) = 1.33_wp
  hhrejfac(1) = 0.4_wp
  hhaccfac(2,1) = One
  hhaccfac(3,1) = OneHalf
  hhaccfac(4,1) = One
  zzacclim(2,2) = 0.8_wp
  zzacclim(3,2) = 1.25_wp
  hhtol(1) = 1.0e-10_wp
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = 0.4_wp
  hhaccfac(3,2) = 1.2_wp
  hhtol(2) = 1.0e-10_wp
end if
if (imethod == 3) then
  ! 'TRUSTOPT'
  ! << TOLS_CVB common block: >>
  ! Criteria for entering local region:
  grd(1,1) = 5.0e-6_wp
  grd(1,2) = 5.0e-6_wp
  sgn(1) = smaller
  sgn(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  dx(1,3) = 5.0e-6_wp
  dx(1,4) = 1.0e-4_wp
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .true.
  nopth1(1) = 5
  nopth2(1) = 2
  delopth1(1) = 0.33333_wp
  delopth2(1) = One
  hhrejfac(1) = 0.08333_wp
  hhaccfac(3,1) = One
  hhtol(1) = 5.0e-6_wp
  dfxmin(1) = Zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = Half
  hhaccfac(3,2) = 1.2_wp
  hhtol(2) = 5.0e-6_wp
  dfxmin(2) = Zero
end if
if (imethod == 4) then
  ! 'DAVIDSON'
  resthr = 1.0e-6_wp
end if
if (imethod == 5) then
  ! 'STEEP'
  exp12tol = Zero
  ! Criteria for entering local region:
  grd(1,1) = 5.0e-6_wp
  grd(1,2) = 5.0e-6_wp
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  ! << TRST_CVB common block: >>
  hhstart = 0.1_wp
  scalesmall(1) = .true.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = Zero
  zzrejmax(1) = 1.33_wp
  hhrejfac(1) = Half
  hhaccfac(2,1) = 1.2_wp
  hhaccfac(3,1) = 1.5_wp
  hhaccfac(4,1) = 1.2_wp
  zzacclim(2,1) = 0.8_wp
  zzacclim(3,1) = 1.25_wp
  hhtol(1) = 5.0e-6_wp
  hhmax(1) = 0.1_wp
  scalesmall(2) = .true.
  nopth1(2) = 1
  nopth2(2) = 0
  zzrejmin(2) = Zero
  zzrejmax(2) = 1.33_wp
  hhrejfac(2) = Half
  hhaccfac(2,2) = 1.2_wp
  hhaccfac(3,2) = 1.5_wp
  hhaccfac(4,2) = 1.2_wp
  zzacclim(2,2) = 0.8_wp
  zzacclim(3,2) = 1.25_wp
  hhtol(2) = 5.0e-6_wp
  hhmax(2) = Half
end if
if ((imethod == 6) .or. (imethod == 7) .or. (imethod == 8) .or. (imethod == 10) .or. (imethod == 12)) then
  ! 'VB2CAS' or 'AUGHESS' or 'AUG2' or 'DFLETCH' or 'SUPER'
  exp12tol = Zero
  ! Criteria for entering local region:
  grd(1,1) = 5.0e-4_wp
  grd(1,2) = 5.0e-4_wp
  sgn(1) = smaller
  sgn(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  grd(1,3) = 5.0e-6_wp
  grd(1,4) = 5.0e-6_wp
  !vv
  !dx(1,3) = 5.0e-6_wp
  dx(1,3) = 5.0e-5_wp
  dx(1,4) = 1.0e-4_wp
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = Zero
  hhrejfac(1) = 0.4_wp
  hhaccfac(2,1) = One
  hhaccfac(3,1) = OneHalf
  hhaccfac(4,1) = One
  zzacclim(2,2) = 0.8_wp
  zzacclim(3,2) = 1.25_wp
  hhtol(1) = 1.0e-10_wp
  dfxmin(1) = Zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = 0.4_wp
  hhaccfac(3,2) = 1.2_wp
  hhtol(2) = 1.0e-10_wp
  dfxmin(2) = Zero
end if

return

end subroutine tunedefs2_cvb
