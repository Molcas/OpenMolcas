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

subroutine initfiles(length,lenv,lenn)
! this routine distributes work space WRK for required files
! for fix mediates it defines also mapd and mapi, for help mediates
! it estimates their length and distributes WRK (i.e. def poss0 parameters)
!
! length - total length of all work space needed
! lenv   - length of V - type array
! lenn   - length of N - type array
!
! !N.B. This routine cannot run with +OP2 level

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: length, lenv, lenn
#include "ccsd1.fh"
#include "ccsd2.fh"
integer(kind=iwp) :: lengthh, lengthm, lengthn, lengthv, maxnoa, maxnorb, maxnvb, maxov(8), posst, symp, sympq, symq, symr, syms

!1 maps and positions for fix mediated

!1.0 maps for DP - diagonal part
!    N.B. DP has one degree of freedom, while other 1 index has none
!    DP1 - dp(p)a
!    DP2 - dp(p)b

do symp=1,nsym
  do symq=1,nsym
    do symr=1,nsym
      mapidp1(symp,symq,symr) = 0
      mapidp2(symp,symq,symr) = 0
    end do
  end do
end do

posst = 1

possdp10 = posst
mapddp1(0,1) = 5
mapddp1(0,2) = 0
mapddp1(0,3) = 0
mapddp1(0,4) = 0
mapddp1(0,5) = nsym
mapddp1(0,6) = 0

do symp=1,nsym
  mapddp1(symp,1) = posst
  mapddp1(symp,2) = norb(symp)
  mapddp1(symp,3) = symp
  mapddp1(symp,4) = 1
  mapddp1(symp,5) = 1
  mapddp1(symp,6) = 1
  mapidp1(symp,1,1) = symp
  posst = posst+norb(symp)
end do

possdp20 = posst
mapddp2(0,1) = 5
mapddp2(0,2) = 0
mapddp2(0,3) = 0
mapddp2(0,4) = 0
mapddp2(0,5) = nsym
mapddp2(0,6) = 0

do symp=1,nsym
  mapddp2(symp,1) = posst
  mapddp2(symp,2) = norb(symp)
  mapddp2(symp,3) = symp
  mapddp2(symp,4) = 1
  mapddp2(symp,5) = 1
  mapddp2(symp,6) = 1
  mapidp2(symp,1,1) = symp
  posst = posst+norb(symp)
end do

!1.1 maps for T1
!    T11 - t1oaa(a,i)
!    T12 - t1obb(a,i)
!    T13 - t1naa(a,i)
!    T14 - t1nbb(a,i)

posst110 = posst
call grc0(2,0,3,1,0,0,1,posst110,posst,mapdt11,mapit11)
posst120 = posst
call grc0(2,0,4,2,0,0,1,posst120,posst,mapdt12,mapit12)
posst130 = posst
call grc0(2,0,3,1,0,0,1,posst130,posst,mapdt13,mapit13)
posst140 = posst
call grc0(2,0,4,2,0,0,1,posst140,posst,mapdt14,mapit14)

!1.2 maps for F1
!    F11 - FI(a,e)aa
!    F12 - FI(a,e)bb

possf110 = posst
call grc0(2,0,3,3,0,0,1,possf110,posst,mapdf11,mapif11)
possf120 = posst
call grc0(2,0,4,4,0,0,1,possf120,posst,mapdf12,mapif12)

!1.3 maps for F2
!    F21 - FII(m,i)aa
!    F22 - FII(m,i)bb

possf210 = posst
call grc0(2,0,1,1,0,0,1,possf210,posst,mapdf21,mapif21)
possf220 = posst
call grc0(2,0,2,2,0,0,1,possf220,posst,mapdf22,mapif22)

!1.4 maps for F3
!    F31 - FIII(e,m)aa
!    F32 - FIII(e,m)bb

possf310 = posst
call grc0(2,0,3,1,0,0,1,possf310,posst,mapdf31,mapif31)
possf320 = posst
call grc0(2,0,4,2,0,0,1,possf320,posst,mapdf32,mapif32)

!1.5 maps for FK
!    FK1 - f(a,b)aa
!    FK2 - f(a,b)bb
!    FK3 - f(a,i)aa
!    FK4 - f(a,i)bb
!    FK5 - f(i,j)aa
!    FK6 - f(i,j)bb

possfk10 = posst
call grc0(2,0,3,3,0,0,1,possfk10,posst,mapdfk1,mapifk1)
possfk20 = posst
call grc0(2,0,4,4,0,0,1,possfk20,posst,mapdfk2,mapifk2)
possfk30 = posst
call grc0(2,0,3,1,0,0,1,possfk30,posst,mapdfk3,mapifk3)
possfk40 = posst
call grc0(2,0,4,2,0,0,1,possfk40,posst,mapdfk4,mapifk4)
possfk50 = posst
call grc0(2,0,1,1,0,0,1,possfk50,posst,mapdfk5,mapifk5)
possfk60 = posst
call grc0(2,0,2,2,0,0,1,possfk60,posst,mapdfk6,mapifk6)

!1.6 maps for T2
!    T21 - t2n(ab,ij)aaaa
!    T22 - t2n(ab,ij)bbbb
!    T33 - t2n(a,b,i,j)abab

posst210 = posst
call grc0(4,4,3,3,1,1,1,posst210,posst,mapdt21,mapit21)
posst220 = posst
call grc0(4,4,4,4,2,2,1,posst220,posst,mapdt22,mapit22)
posst230 = posst
call grc0(4,0,3,4,1,2,1,posst230,posst,mapdt23,mapit23)

!1.7 maps for W0
!    W01 - <mn||ij>aaaa
!    W02 - <mn||ij>bbbb
!    W03 - <mn||ij>abab

possw010 = posst
call grc0(4,4,1,1,1,1,1,possw010,posst,mapdw01,mapiw01)
possw020 = posst
call grc0(4,4,2,2,2,2,1,possw020,posst,mapdw02,mapiw02)
possw030 = posst
call grc0(4,0,1,2,1,2,1,possw030,posst,mapdw03,mapiw03)

!1.8 maps for W1
!    W11 - <ie||mn>aaaa
!    W12 - <ie||mn>bbbb
!    W13 - <ie||mn>abab
!    W14 - <ie||mn>baab

possw110 = posst
call grc0(4,3,1,3,1,1,1,possw110,posst,mapdw11,mapiw11)
possw120 = posst
call grc0(4,3,2,4,2,2,1,possw120,posst,mapdw12,mapiw12)
possw130 = posst
call grc0(4,0,1,4,1,2,1,possw130,posst,mapdw13,mapiw13)
possw140 = posst
call grc0(4,0,2,3,1,2,1,possw140,posst,mapdw14,mapiw14)

!2 for help files mapps are irrelevant,
!  here only estimation of maximal length is done to
!  define poss0 of help files
!  we have:
!  four V files - of vvoo type
!  four M files - of vvo  type
!  four H files - of voo  type
!  one  N file  - of nn   type

!2.* def max{noa}, max{norb}, max{nvb}, maxov(isym)=max{noa(isym),nvb(isym)}

maxnoa = noa(1)
maxnvb = nvb(1)
maxnorb = norb(1)
do symp=1,nsym
  if (noa(symp) > maxnoa) maxnoa = noa(symp)
  if (norb(symp) > maxnorb) maxnorb = norb(symp)
  if (nvb(symp) > maxnvb) maxnvb = nvb(symp)
  if (nvb(symp) > noa(symp)) then
    maxov(symp) = nvb(symp)
  else
    maxov(symp) = noa(symp)
  end if
end do

!2.* def lengths of V,M,H and N files

lengthv = 0
lengthm = 0
lengthh = 0
lengthn = 0

do symp=1,nsym
  ! symq is not known for N file
  ! instead of norb(symr) maxnorb will be used so that reallength<=length
  lengthn = lengthn+norb(symp)*maxnorb
  do symq=1,nsym
    sympq = mmul(symp,symq)
    ! symr is not known for M and H files
    ! instead of noa(symr) maxnoa will be used so that reallength<=length
    lengthm = lengthm+maxov(symp)*maxov(symq)*maxnoa
    lengthh = lengthh+maxov(symp)*noa(symq)*maxnoa
    do symr=1,nsym
      syms = mmul(sympq,symr)
      lengthv = lengthv+maxov(symp)*maxov(symq)*noa(symr)*noa(syms)
    end do
  end do
end do

!2.1 V - files

possv10 = posst
posst = posst+lengthv
possv20 = posst
posst = posst+lengthv
possv30 = posst
posst = posst+lengthv
possv40 = posst
posst = posst+lengthv
lenv = lengthv

!2.2 M - files

possm10 = posst
posst = posst+lengthm
possm20 = posst
posst = posst+lengthm
possm30 = posst
posst = posst+lengthm
possm40 = posst
posst = posst+lengthm

!2.3 H - files

possh10 = posst
posst = posst+lengthh
possh20 = posst
posst = posst+lengthh
possh30 = posst
posst = posst+lengthh
possh40 = posst
posst = posst+lengthh

!2.4 N,P - files

possn0 = posst
posst = posst+lengthn
possp0 = posst
posst = posst+lengthn
lenn = lengthn

!2.5 declare space for help matrix D in for matrix multiplication C=AT*B if
!    mchntyp=2

if (mchntyp == 2) then
  possd0 = posst
  if (maxnoa <= maxnvb) then
    posst = posst+maxnoa*maxnoa*maxnvb*maxnvb
  else
    posst = posst+maxnoa*maxnoa*maxnoa*maxnoa
  end if
end if

!2.6 def size of Work space
length = posst-1

return

end subroutine initfiles
