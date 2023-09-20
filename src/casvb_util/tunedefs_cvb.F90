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

implicit real*8(a-h,o-z)
save huge, smallest
data huge/1d20/
data smallest/1d-10/

! General defaults
! << TUNE_CVB common block: >>
cnrmtol = 1.d-9
safety = 1d-7
signtol = 1.d-3
alftol = 1.d-10
dfxtol = 1.d-10
exp12tol = huge
grdwrngtol = -huge
eigwrngtol = -huge
endwhenclose = .false.
! << TOLS_CVB common block: >>
! (*,1) ... global region, non-singular Hessian
! (*,2) ... global region, singular Hessian
! (*,3) ... local region, non-singular Hessian
! (*,4) ... local region, singular Hessian
! (*,5) ... wrong stationary point, non-singular Hessian
! (*,6) ... wrong stationary point, singular Hessian
! First set values that disable tests:
singul(1) = -huge
singul(2) = -huge
singul(3) = -huge
do j=1,6
  dfx(j) = huge
  sgn(j) = huge
  zzmax(j) = huge
  zzmin(j) = -huge
  do i=1,3
    dx(i,j) = huge
    grd(i,j) = huge
  end do
end do
! << TRST_CVB common block: >>
scalesmall(1) = .false.
nopth1(1) = 1
nopth2(1) = 0
delopth1(1) = 1d0
delopth2(1) = huge
delopth1(2) = 1d0
delopth2(2) = huge
hhrejfac(1) = 1d0
hhaccfac(1,1) = 1d0
hhaccfac(2,1) = 1d0
hhaccfac(3,1) = 1d0
hhaccfac(4,1) = 1d0
hhaccfac(5,1) = 1d0
zzacclim(1,1) = -huge
zzacclim(2,1) = -huge
zzacclim(3,1) = huge
zzacclim(4,1) = huge
hhtol(1) = -huge
hhmax(1) = 1d0
dfxmin(1) = -huge
zzrejmin(1) = -huge
zzrejmax(1) = huge
scalesmall(2) = .false.
nopth1(2) = 1
nopth2(2) = 0
hhrejfac(2) = 1d0
hhaccfac(1,2) = 1d0
hhaccfac(2,2) = 1d0
hhaccfac(3,2) = 1d0
hhaccfac(4,2) = 1d0
hhaccfac(5,2) = 1d0
zzacclim(1,2) = -huge
zzacclim(2,2) = -huge
zzacclim(3,2) = huge
zzacclim(4,2) = huge
hhtol(2) = -huge
hhmax(2) = 1d0
dfxmin(2) = -huge
zzrejmin(2) = -huge
zzrejmax(2) = huge
hhstart = 1d0
! << DAVTUNE global vars: >>
resthr = 5d-6
orththr = smallest
nortiter = 50
mxdav = 200
follow = .false.

return

end subroutine tunedefs_cvb
