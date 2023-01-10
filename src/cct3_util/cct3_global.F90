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

module CCT3_global

! Common variables of whole T3 program
!
!1    max. number of orbitals (basis functions)
!     - maxorb
!
!2    characteristics from MOLCAS
!
!2.1  Number of active electrons
!     - nactel (*)
!2.2  spin state of the system
!     - ispin
!2.3  number of irreps in the system
!     - nsym
!2.4  symmetry state of the system
!     - lsym
!2.5  matrix multiplication table
!     - mmul
!2.6  vectors containing size of the system
!     - noa,nob,nva,nvb,norb
!2.7  orbital energies
!     - eps
!
!3    internal CCSD characteristics
!
!3.1  size characteristic Table
!     - dimm
!3.2  shift vector
!     - nshf
!
!4    input parameters (read from input file)
!
!4.1  title of the job + number ot tilte lines
!     - title,ntit (*)
!4.3  type of t3 contribution
!     - typt3
!4.4  type of denominator
!     - typden
!4.10 spin adaptation key
!     - keysa
!4.12 restart information file name
!     - filerst
!4.13 type of machine (from point of mtx multiplication efficiency)
!     (parameter posd0 will be defined in initfile if mchntyp=2)
!     - mchntyp,posd0
!4.14 limiting ratio for using AT*B mtx multiplication, if mchntyp=2
!     - slim
!4.15 denominator shifts
!     - shifto,shiftv
!4.16 maximal allowed work space
!     - maxspace
!4.17 level of printing control parameter
!     - fullprint
!4.17 No operation key
!     - noop
!4.18 I/O control key
!     - iokey
!4.19 Matrix handling control key
!     - mhkey
!4.20 Key for swithching segmentation of I,J cycle
!     - ijsegkey
!4.21 Lower and upper limitations of I,J cycle (for ijsegkey=1)
!     - symimin,symimax,symjmin,symjmax
!     - imin,imax,jmin,jmax
!
!5    diska address file
!     - daddr
!
! Maps for all mediates, used in NIT3
!
!1   maps for fixed mediates
!
!1.0 maps for DP - diagonal part
!     DP1 - dp(p)a
!     DP2 - dp(p)b
!
!1.1 maps for T1
!     T11 - t1oaa(a,i)
!     T12 - t1obb(a,i)
!
!1.5 maps for FK
!     FK1 - f(a,b)aa
!     FK2 - f(a,b)bb
!     FK3 - f(a,i)aa
!     FK4 - f(a,i)bb
!     FK5 - f(i,j)aa
!     FK6 - f(i,j)bb
!
!1.6 maps for T2
!     T21 - t2n(ab,ij)aaaa
!     T22 - t2n(ab,ij)bbbb
!     T33 - t2n(a,b,i,j)abab
!
!1.8 maps for W1
!     W11 - <ie||mn>aaaa
!     W12 - <ie||mn>bbbb
!     W13 - <ie||mn>abab
!     W14 - <ie||nm>baab
!
!1.9 maps for W2
!     W21 - <ab||ij>aaaa
!     W22 - <ab||ij>bbbb
!     W23 - <a,b|i,j>abab
!
!2   maps for help files
!    There are :
!     2  WX,VX files - vv2 type
!     2  L     files - vvv(vvo) type
!     3  RX    files - vv2* type
!     3  M     files - of vv (vo)  type
!     3  H     files - of v (o)  type
!     2  N,PX  files - of nn   type

use Definitions, only: wp, iwp

implicit none
private

! pos0 - initial position
! d    - direct map
! i    - inverse map
type Map_Type
  integer(kind=iwp) :: d(0:512,6), i(8,8,8), pos0
end type Map_Type

integer(kind=iwp), parameter :: maxorb = 1024
character(len=*), parameter :: t3nam = 'T3VVVO'
integer(kind=iwp) :: daddr(8), dimm(5,8), fullprint, ijsegkey, imax, imin, iokey, ispin, jmax, jmin, keysa, lsym, maxspace, &
                     mchntyp, mhkey, mmul(8,8), noa(8), nob(8), noop, norb(8), nshf(maxorb), nsym, nva(8), nvb(8), posd0, &
                     symimax, symimin, symjmax, symjmin, T3IntPos(maxorb), T3Off(512,8), typden, typt3
real(kind=wp) :: eps(maxorb), shifto, shiftv, slim
character(len=6) :: filerst
type(Map_Type) :: dp1, dp2, fk1, fk2, fk3, fk4, fk5, fk6, h1, h2, h3, l1, l2, m1, m2, m3, n, px, rx1, rx2, rx3, t11, t12, t21, &
                  t22, t23, vx, w11, w12, w13, w14, w21, w22, w23, wx

public :: daddr, dimm, dp1, dp2, eps, filerst, fk1, fk2, fk3, fk4, fk5, fk6, fullprint, h1, h2, h3, ijsegkey, imax, imin, iokey, &
          ispin, jmax, jmin, keysa, l1, l2, lsym, m1, m2, m3, Map_Type, maxorb, maxspace, mchntyp, mhkey, mmul, n, noa, nob, noop, &
          norb, nshf, nsym, nva, nvb, posd0, px, rx1, rx2, rx3, shifto, shiftv, slim, symimax, symimin, symjmax, symjmin, t11, &
          t12, t21, t22, t23, typden, typt3, vx, w11, w12, w13, w14, w21, w22, w23, wx

public :: T3IntPos, t3nam, T3Off

end module CCT3_global
