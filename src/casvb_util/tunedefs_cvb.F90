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

subroutine tunedefs_cvb()

use casvb_global, only: alftol, cnrmtol, delopth1, delopth2, dfx, dfxmin, dfxtol, dx, eigwrngtol, endwhenclose, exp12tol, follow, &
                        grd, grdwrngtol, hhaccfac, hhmax, hhrejfac, hhstart, hhtol, mxdav, nopth1, nopth2, nortiter, orththr, &
                        resthr, safety, scalesmall, sgn, signtol, singul, zzacclim, zzmax, zzmin, zzrejmax, zzrejmin

use Constants, only: One
use Definitions, only: wp

implicit none
real(kind=wp), parameter :: hge = 1.0e20_wp, smallest = 1.0e-10_wp

! General defaults
! << TUNE_CVB common block: >>
cnrmtol = 1.0e-9_wp
safety = 1.0e-7_wp
signtol = 1.0e-3_wp
alftol = 1.0e-10_wp
dfxtol = 1.0e-10_wp
exp12tol = hge
grdwrngtol = -hge
eigwrngtol = -hge
endwhenclose = .false.
! << TOLS_CVB common block: >>
! (*,1) ... global region, non-singular Hessian
! (*,2) ... global region, singular Hessian
! (*,3) ... local region, non-singular Hessian
! (*,4) ... local region, singular Hessian
! (*,5) ... wrong stationary point, non-singular Hessian
! (*,6) ... wrong stationary point, singular Hessian
! First set values that disable tests:
singul(1) = -hge
singul(2) = -hge
singul(3) = -hge
dfx(:) = hge
sgn(:) = hge
zzmax(:) = hge
zzmin(:) = -hge
dx(:,:) = hge
grd(:,:) = hge
! << TRST_CVB common block: >>
scalesmall(1) = .false.
nopth1(1) = 1
nopth2(1) = 0
delopth1(1) = One
delopth2(1) = hge
delopth1(2) = One
delopth2(2) = hge
hhrejfac(1) = One
hhaccfac(1,1) = One
hhaccfac(2,1) = One
hhaccfac(3,1) = One
hhaccfac(4,1) = One
hhaccfac(5,1) = One
zzacclim(1,1) = -hge
zzacclim(2,1) = -hge
zzacclim(3,1) = hge
zzacclim(4,1) = hge
hhtol(1) = -hge
hhmax(1) = One
dfxmin(1) = -hge
zzrejmin(1) = -hge
zzrejmax(1) = hge
scalesmall(2) = .false.
nopth1(2) = 1
nopth2(2) = 0
hhrejfac(2) = One
hhaccfac(1,2) = One
hhaccfac(2,2) = One
hhaccfac(3,2) = One
hhaccfac(4,2) = One
hhaccfac(5,2) = One
zzacclim(1,2) = -hge
zzacclim(2,2) = -hge
zzacclim(3,2) = hge
zzacclim(4,2) = hge
hhtol(2) = -hge
hhmax(2) = One
dfxmin(2) = -hge
zzrejmin(2) = -hge
zzrejmax(2) = hge
hhstart = One
! << DAVTUNE global vars: >>
resthr = 5.0e-6_wp
orththr = smallest
nortiter = 50
mxdav = 200
follow = .false.

return

end subroutine tunedefs_cvb
