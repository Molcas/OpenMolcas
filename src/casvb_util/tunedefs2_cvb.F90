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

subroutine tunedefs2_cvb(imethod,endifclose1)

implicit real*8(a-h,o-z)
logical endifclose1
#include "tune_cvb.fh"
#include "tols_cvb.fh"
#include "trst_cvb.fh"
#include "davtune_cvb.fh"
save small, small2, smaller, zero
data small/1d-3/,small2/1d-5/,smaller/1d-6/
data zero/0d0/

! << TUNE_CVB common block: >>
endifclose = endifclose1
if ((imethod == 1) .or. (imethod == 10)) then
  !'FLETCHER'
  exp12tol = zero
  ! Criteria for entering local region:
  grd(1,1) = 5d-4
  grd(1,2) = 5d-4
  sign(1) = smaller
  sign(2) = smaller
  singul(1) = 1d-2
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  grd(1,3) = 5d-6
  grd(1,4) = 5d-6
  dx(1,3) = 5d-6
  dx(1,4) = 1d-4
  singul(2) = 1d-2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = 0d0
  hhrejfac(1) = .4d0
  hhaccfac(2,1) = 1d0
  hhaccfac(3,1) = 1.5d0
  hhaccfac(4,1) = 1d0
  zzacclim(2,2) = .8d0
  zzacclim(3,2) = 1.25d0
  hhtol(1) = 1d-10
  dfxmin(1) = zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = .4d0
  hhaccfac(3,2) = 1.2d0
  hhtol(2) = 1d-10
  dfxmin(2) = zero
end if
if (imethod == 2) then
  ! 'TRIM'
  exp12tol = zero
  ! Criteria for entering local region:
  grd(1,1) = 5d-6
  grd(1,2) = 5d-6
  sign(1) = smaller
  sign(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  dx(1,3) = 5d-6
  dx(1,4) = 1d-4
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = .75d0
  zzrejmax(1) = 1.33d0
  hhrejfac(1) = .4d0
  hhaccfac(2,1) = 1d0
  hhaccfac(3,1) = 1.5d0
  hhaccfac(4,1) = 1d0
  zzacclim(2,2) = .8d0
  zzacclim(3,2) = 1.25d0
  hhtol(1) = 1d-10
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = .4d0
  hhaccfac(3,2) = 1.2d0
  hhtol(2) = 1d-10
end if
if (imethod == 3) then
  ! 'TRUSTOPT'
  ! << TOLS_CVB common block: >>
  ! Criteria for entering local region:
  grd(1,1) = 5d-6
  grd(1,2) = 5d-6
  sign(1) = smaller
  sign(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  dx(1,3) = 5d-6
  dx(1,4) = 1d-4
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .true.
  nopth1(1) = 5
  nopth2(1) = 2
  delopth1(1) = .33333d0
  delopth2(1) = 1.d0
  hhrejfac(1) = .08333d0
  hhaccfac(3,1) = 1d0
  hhtol(1) = 5.d-6
  dfxmin(1) = zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = .5d0
  hhaccfac(3,2) = 1.2d0
  hhtol(2) = 5.d-6
  dfxmin(2) = zero
end if
if (imethod == 4) then
  ! 'DAVIDSON'
  resthr = 1d-6
end if
if (imethod == 5) then
  ! 'STEEP'
  exp12tol = zero
  ! Criteria for entering local region:
  grd(1,1) = 5d-6
  grd(1,2) = 5d-6
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  ! << TRST_CVB common block: >>
  hhstart = .1d0
  scalesmall(1) = .true.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = 0d0
  zzrejmax(1) = 1.33d0
  hhrejfac(1) = .5d0
  hhaccfac(2,1) = 1.2d0
  hhaccfac(3,1) = 1.5d0
  hhaccfac(4,1) = 1.2d0
  zzacclim(2,1) = .8d0
  zzacclim(3,1) = 1.25d0
  hhtol(1) = 5.d-6
  hhmax(1) = .1d0
  scalesmall(2) = .true.
  nopth1(2) = 1
  nopth2(2) = 0
  zzrejmin(2) = 0d0
  zzrejmax(2) = 1.33d0
  hhrejfac(2) = .5d0
  hhaccfac(2,2) = 1.2d0
  hhaccfac(3,2) = 1.5d0
  hhaccfac(4,2) = 1.2d0
  zzacclim(2,2) = .8d0
  zzacclim(3,2) = 1.25d0
  hhtol(2) = 5.d-6
  hhmax(2) = .5d0
end if
if ((imethod == 6) .or. (imethod == 7) .or. (imethod == 8) .or. (imethod == 10) .or. (imethod == 12)) then
  ! 'VB2CAS' or 'AUGHESS' or 'AUG2' or 'DFLETCH' or 'SUPER'
  exp12tol = zero
  ! Criteria for entering local region:
  grd(1,1) = 5d-4
  grd(1,2) = 5d-4
  sign(1) = smaller
  sign(2) = smaller
  singul(1) = small
  zzmin(1) = -small
  zzmin(2) = -small
  ! Final convergence criteria:
  grd(1,3) = 5d-6
  grd(1,4) = 5d-6
  !vv
  !dx(1,3) = 5d-6
  dx(1,3) = 5d-5
  dx(1,4) = 1d-4
  singul(2) = small2
  ! << TRST_CVB common block: >>
  scalesmall(1) = .false.
  nopth1(1) = 1
  nopth2(1) = 0
  zzrejmin(1) = 0d0
  hhrejfac(1) = .4d0
  hhaccfac(2,1) = 1d0
  hhaccfac(3,1) = 1.5d0
  hhaccfac(4,1) = 1d0
  zzacclim(2,2) = .8d0
  zzacclim(3,2) = 1.25d0
  hhtol(1) = 1d-10
  dfxmin(1) = zero
  scalesmall(2) = .false.
  nopth1(2) = 1
  nopth2(2) = 0
  hhrejfac(2) = .4d0
  hhaccfac(3,2) = 1.2d0
  hhtol(2) = 1d-10
  dfxmin(2) = zero
end if

return

end subroutine tunedefs2_cvb
