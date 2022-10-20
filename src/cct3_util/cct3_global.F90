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
!     - nactel
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
!     - title,ntit
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
!     2  W,V  files - vv2 type
!     2    L  files - vvv(vvo) type
!     3    R  files - vv2* type
!     3    M  files - of vv (vo)  type
!     3    H  files - of v (o)  type
!     2   N,P files - of nn   type

use Definitions, only: wp, iwp

implicit none
private

integer(kind=iwp), parameter :: maxorb = 1024
integer(kind=iwp) :: daddr(8), dimm(5,8), fullprint, ijsegkey, imax, imin, iokey, ispin, jmax, jmin, keysa, lsym, maxspace, &
                     mchntyp, mhkey, mmul(8,8), nactel, noa(8), nob(8), noop, norb(8), nshf(maxorb), nsym, ntit, nva(8), nvb(8), &
                     posd0, symimax, symimin, symjmax, symjmin, typden, typt3
real(kind=wp) :: eps(maxorb), shifto, shiftv, slim
character(len=72) :: title
character(len=6) :: filerst

public :: daddr, dimm, eps, filerst, fullprint, ijsegkey, imax, imin, iokey, ispin, jmax, jmin, keysa, lsym, maxorb, maxspace, &
          mchntyp, mhkey, mmul, nactel, noa, nob, noop, norb, nshf, nsym, ntit, nva, nvb, posd0, shifto, shiftv, slim, symimax, &
          symimin, symjmax, symjmin, title, typden, typt3

#define _MAP_D_ 0:512,6
#define _MAPI_ 8,8,8
integer(kind=iwp), public :: mapddp1(_MAP_D_), mapddp2(_MAP_D_), mapdfk1(_MAP_D_), mapdfk2(_MAP_D_), mapdfk3(_MAP_D_), &
                             mapdfk4(_MAP_D_), mapdfk5(_MAP_D_), mapdfk6(_MAP_D_), mapdh1(_MAP_D_), mapdh2(_MAP_D_), &
                             mapdh3(_MAP_D_), mapdl1(_MAP_D_), mapdl2(_MAP_D_), mapdm1(_MAP_D_), mapdm2(_MAP_D_), mapdm3(_MAP_D_), &
                             mapdn(_MAP_D_), mapdp(_MAP_D_), mapdr1(_MAP_D_), mapdr2(_MAP_D_), mapdr3(_MAP_D_), mapdt11(_MAP_D_), &
                             mapdt12(_MAP_D_), mapdt21(_MAP_D_), mapdt22(_MAP_D_), mapdt23(_MAP_D_), mapdv(_MAP_D_), &
                             mapdw(_MAP_D_), mapdw11(_MAP_D_), mapdw12(_MAP_D_), mapdw13(_MAP_D_), mapdw14(_MAP_D_), &
                             mapdw21(_MAP_D_), mapdw22(_MAP_D_), mapdw23(_MAP_D_), &
                             mapidp1(_MAPI_), mapidp2(_MAPI_), mapifk1(_MAPI_), mapifk2(_MAPI_), mapifk3(_MAPI_), mapifk4(_MAPI_), &
                             mapifk5(_MAPI_), mapifk6(_MAPI_), mapih1(_MAPI_), mapih2(_MAPI_), mapih3(_MAPI_), mapil1(_MAPI_), &
                             mapil2(_MAPI_), mapim1(_MAPI_), mapim2(_MAPI_), mapim3(_MAPI_), mapin(_MAPI_), mapip(_MAPI_), &
                             mapir1(_MAPI_), mapir2(_MAPI_), mapir3(_MAPI_), mapit11(_MAPI_), mapit12(_MAPI_), mapit21(_MAPI_), &
                             mapit22(_MAPI_), mapit23(_MAPI_), mapiv(_MAPI_), mapiw(_MAPI_), mapiw11(_MAPI_), mapiw12(_MAPI_), &
                             mapiw13(_MAPI_), mapiw14(_MAPI_), mapiw21(_MAPI_), mapiw22(_MAPI_), mapiw23(_MAPI_), &
                             posdp10, posdp20, posfk10, posfk20, posfk30, posfk40, posfk50, posfk60, posh10, posh20, posh30, &
                             posl10, posl20, posm10, posm20, posm30, posn0, posp0, posr10, posr20, posr30, post110, post120, &
                             post210, post220, post230, posv0, posw0, posw110, posw120, posw130, posw140, posw210, posw220, posw230

end module CCT3_global
