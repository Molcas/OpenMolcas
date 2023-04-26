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

module ccsd_global

! A   Most common arrays of whole program
!
!1    max. number of orbitals (basis functions)
!     - maxorb
!
!2    characteristics from MOLCAS
!
!2.1  Number of active electrons
!     - nactel
!
!2.2  spin state of the system
!     - ispin
!
!2.3  number of irreps in the system
!     - nsym
!
!2.4  symmetry state of the system
!     - lsym
!
!2.5  matrix multiplication table
!     - mmul
!
!2.6  vectors containing size of the system
!     - noa, nob, nva, nvb, norb
!
!2.7  orbital energies
!     - eps
!
!3    internal CCSD characteristics
!
!3.1  size characteristic Table
!     - dimm(1:5,1:8)
!
!3.2  shift vector
!     - nshf(1:maxorb)
!
!4    input parameters (read from input file)
!
!4.1  title of the job + number ot title lines
!     - title, ntit
!
!4.2  Maximal number of iterations
!     - maxiter
!
!4.3  type of t3 contribution
!     - typt3
!
!4.4  type of denominator
!     - typden
!
!4.5  using of extrapolation
!     - yesext
!
!4.6  first iteration of extrapolation
!     - firstext
!
!4.7  extrapolation cycle
!     - cycext
!
!4.9  tolerance for energy
!     - ccconv
!
!4.10 spin adaptation key
!     - keysa
!
!4.11 restart key
!     - keyrst
!
!4.12 restart information file name
!     - filerst
!
!4.13 type of machine (from point of mtx multiplication efficiency)
!     parameter posd0 will be defined in initfile if mchntyp=2)
!     - mchntyp, posd0
!
!4.14 limiting ratio for using AT*B mtx multiplication, if mchntyp=2
!     - slim
!
!4.15 denominator shifts
!     - shifto, shiftv
!
!4.16 maximal allowed work space
!     - maxspace
!
!4.17 level of printing option
!     - fullprint
!
!4.18 No Operation key
!     - noop
!
!4.19 Input/Output style key
!     - iokey
!
!4.20 Matrix handling control key
!     - mhkey
!
!4.21 Key to switch off CCSD part
!     - noccsd
!
!FUE  added to transport the total energies
!     - Escf
!
! B   Maps for all mediates, used in CC
!
!1    maps for fixed mediates
!
!1.0  maps for DP - diagonal part
!     DP1 - dp(p)a
!     DP2 - dp(p)b
!
!1.1  maps for T1
!     T11 - t1oaa(a,i)
!     T12 - t1obb(a,i)
!     T13 - t1naa(a,i)
!     T14 - t1nbb(a,i)
!
!1.2  maps for F1
!     F11 - FI(a,e)aa
!     F12 - FI(a,e)bb
!
!1.3  maps for F2
!     F21 - FII(m,i)aa
!     F22 - FII(m,i)bb
!
!1.4  maps for F3
!     F31 - FIII(e,m)aa
!     F32 - FIII(e,m)bb
!
!1.5  maps for FK
!     FK1 - f(a,b)aa
!     FK2 - f(a,b)bb
!     FK3 - f(a,i)aa
!     FK4 - f(a,i)bb
!     FK5 - f(i,j)aa
!     FK6 - f(i,j)bb
!
!1.6  maps for T2
!     T21 - t2n(ab,ij)aaaa
!     T22 - t2n(ab,ij)bbbb
!     T33 - t2n(a,b,i,j)abab
!
!1.7  maps for W0
!     W01 - <mn||ij>aaaa
!     W02 - <mn||ij>bbbb
!     W03 - <mn||ij>abab
!
!1.8  maps for W1
!     W11 - <ie||mn>aaaa
!     W12 - <ie||mn>bbbb
!     W13 - <ie||mn>abab
!     W14 - <ie||nm>baab
!
!2    maps for help files
!
!     There are :
!     four V  files - of vvoo type
!     four M  files - of vvo  type
!     four H  files - of voo  type
!     two N,P files - of nn   type
!     For V1-V4, M1-M4 H1-H4 and N,P there is corresponding space on WRK at
!     positions %pos0,
!     that are ALWAYS used when using corresponding help files
!     There are also corresponding %d and %i matrix, that are USUALLY used
!     when using corresponding help files
!
!     Moreover, there are 6 additional R1-R6, that can be used
!     arbitrarily, of course with no corresponding space on WRK
!
! C   Parallel run
!
!1    general parameters
!     - maxproc
!
!2    parameters for sumoverab part
!     - nprocab, idab, ideffab, idtmab, ididle
!
!3    parameters for sumoverb and intermezzo part
!     - idaaaa, idbbbb, idaabb, idabba, idbaab, idbbaa
!
!4    parameters for finale part
!     - idfin

use Definitions, only: wp, iwp

implicit none
private

! pos0 - initial position
! d    - direct map
! i    - inverse map
type Map_Type
  integer(kind=iwp) :: d(0:512,6), i(8,8,8), pos0
end type Map_Type

integer(kind=iwp), parameter :: minfiles = 17, maxfiles = 50, &
                                maxorb = 1024, maxproc = 16

integer(kind=iwp) :: cycext, daddr(10:maxfiles) = 0, dimm(5,8), filestatus(10:maxfiles) = 0, firstext, fullprint, idaaaa, idaabb, &
                     idab(maxproc), idabba, idbaab, idbbaa, idbbbb, idfin, iokey, ispin, keyrst, keysa, lsym, maxiter, maxspace, &
                     mchntyp, mhkey, mmul(8,8), nactel, noa(8), nob(8), noccsd, noop, norb(8), nprocab, nshf(maxorb), nsym, ntit, &
                     nva(8), nvb(8), posd0, typden, typt3, yesext
real(kind=wp) :: ccconv, eps(maxorb), Escf, ideffab(maxproc), ididle(maxproc), idtmab(maxproc), shifto, shiftv, slim
character(len=72) :: title
character(len=6) :: filename(10:maxfiles) = 'Temp', filerst
type(Map_Type) :: dp1, dp2, f11, f12, f21, f22, f31, f32, fk1, fk2, fk3, fk4, fk5, fk6, h1, h2, h3, h4, m1, m2, m3, m4, n, p, r1, &
                  r2, r3, r4, r5, r6, t11, t12, t13, t14, t21, t22, t23, v1, v2, v3, v4, w01, w02, w03, w11, w12, w13, w14

public :: ccconv, cycext, daddr, dimm, dp1, dp2, eps, Escf, f11, f12, f21, f22, f31, f32, filename, filerst, filestatus, firstext, &
          fk1, fk2, fk3, fk4, fk5, fk6, fullprint, h1, h2, h3, h4, idaaaa, idaabb, idab, idabba, idbaab, idbbaa, idbbbb, ideffab, &
          idfin, ididle, idtmab, iokey, ispin, keyrst, keysa, lsym, m1, m2, m3, m4, Map_Type, maxfiles, maxiter, maxorb, maxproc, &
          maxspace, mchntyp, mhkey, minfiles, mmul, n, nactel, noa, nob, noccsd, noop, norb, nprocab, nshf, nsym, ntit, nva, nvb, &
          p, posd0, r1, r2, r3, r4, r5, r6, shifto, shiftv, slim, t11, t12, t13, t14, t21, t22, t23, title, typden, typt3, v1, v2, &
          v3, v4, w01, w02, w03, w11, w12, w13, w14, yesext

end module ccsd_global
