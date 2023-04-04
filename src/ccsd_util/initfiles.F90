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

use ccsd_global, only: dp1, dp2, f11, f12, f21, f22, f31, f32, fk1, fk2, fk3, fk4, fk5, fk6, h1, h2, h3, h4, m1, m2, m3, m4, &
                       mchntyp, mmul, n, noa, norb, nsym, nvb, p, possd0, t11, t12, t13, t14, t21, t22, t23, v1, v2, v3, v4, w01, &
                       w02, w03, w11, w12, w13, w14
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: length, lenv, lenn
integer(kind=iwp) :: lengthh, lengthm, lengthn, lengthv, maxnoa, maxnorb, maxnvb, maxov(8), posst, symp, sympq, symq, symr, syms

!1 maps and positions for fix mediated

!1.0 maps for DP - diagonal part
!    N.B. DP has one degree of freedom, while other 1 index has none
!    DP1 - dp(p)a
!    DP2 - dp(p)b

do symp=1,nsym
  do symq=1,nsym
    do symr=1,nsym
      dp1%i(symp,symq,symr) = 0
      dp2%i(symp,symq,symr) = 0
    end do
  end do
end do

posst = 1

dp1%pos0 = posst
dp1%d(0,1) = 5
dp1%d(0,2) = 0
dp1%d(0,3) = 0
dp1%d(0,4) = 0
dp1%d(0,5) = nsym
dp1%d(0,6) = 0

do symp=1,nsym
  dp1%d(symp,1) = posst
  dp1%d(symp,2) = norb(symp)
  dp1%d(symp,3) = symp
  dp1%d(symp,4) = 1
  dp1%d(symp,5) = 1
  dp1%d(symp,6) = 1
  dp1%i(symp,1,1) = symp
  posst = posst+norb(symp)
end do

dp2%pos0 = posst
dp2%d(0,1) = 5
dp2%d(0,2) = 0
dp2%d(0,3) = 0
dp2%d(0,4) = 0
dp2%d(0,5) = nsym
dp2%d(0,6) = 0

do symp=1,nsym
  dp2%d(symp,1) = posst
  dp2%d(symp,2) = norb(symp)
  dp2%d(symp,3) = symp
  dp2%d(symp,4) = 1
  dp2%d(symp,5) = 1
  dp2%d(symp,6) = 1
  dp2%i(symp,1,1) = symp
  posst = posst+norb(symp)
end do

!1.1 maps for T1
!    T11 - t1oaa(a,i)
!    T12 - t1obb(a,i)
!    T13 - t1naa(a,i)
!    T14 - t1nbb(a,i)

t11%pos0 = posst
call grc0(2,0,3,1,0,0,1,t11%pos0,posst,t11%d,t11%i)
t12%pos0 = posst
call grc0(2,0,4,2,0,0,1,t12%pos0,posst,t12%d,t12%i)
t13%pos0 = posst
call grc0(2,0,3,1,0,0,1,t13%pos0,posst,t13%d,t13%i)
t14%pos0 = posst
call grc0(2,0,4,2,0,0,1,t14%pos0,posst,t14%d,t14%i)

!1.2 maps for F1
!    F11 - FI(a,e)aa
!    F12 - FI(a,e)bb

f11%pos0 = posst
call grc0(2,0,3,3,0,0,1,f11%pos0,posst,f11%d,f11%i)
f12%pos0 = posst
call grc0(2,0,4,4,0,0,1,f12%pos0,posst,f12%d,f12%i)

!1.3 maps for F2
!    F21 - FII(m,i)aa
!    F22 - FII(m,i)bb

f21%pos0 = posst
call grc0(2,0,1,1,0,0,1,f21%pos0,posst,f21%d,f21%i)
f22%pos0 = posst
call grc0(2,0,2,2,0,0,1,f22%pos0,posst,f22%d,f22%i)

!1.4 maps for F3
!    F31 - FIII(e,m)aa
!    F32 - FIII(e,m)bb

f31%pos0 = posst
call grc0(2,0,3,1,0,0,1,f31%pos0,posst,f31%d,f31%i)
f32%pos0 = posst
call grc0(2,0,4,2,0,0,1,f32%pos0,posst,f32%d,f32%i)

!1.5 maps for FK
!    FK1 - f(a,b)aa
!    FK2 - f(a,b)bb
!    FK3 - f(a,i)aa
!    FK4 - f(a,i)bb
!    FK5 - f(i,j)aa
!    FK6 - f(i,j)bb

fk1%pos0 = posst
call grc0(2,0,3,3,0,0,1,fk1%pos0,posst,fk1%d,fk1%i)
fk2%pos0 = posst
call grc0(2,0,4,4,0,0,1,fk2%pos0,posst,fk2%d,fk2%i)
fk3%pos0 = posst
call grc0(2,0,3,1,0,0,1,fk3%pos0,posst,fk3%d,fk3%i)
fk4%pos0 = posst
call grc0(2,0,4,2,0,0,1,fk4%pos0,posst,fk4%d,fk4%i)
fk5%pos0 = posst
call grc0(2,0,1,1,0,0,1,fk5%pos0,posst,fk5%d,fk5%i)
fk6%pos0 = posst
call grc0(2,0,2,2,0,0,1,fk6%pos0,posst,fk6%d,fk6%i)

!1.6 maps for T2
!    T21 - t2n(ab,ij)aaaa
!    T22 - t2n(ab,ij)bbbb
!    T33 - t2n(a,b,i,j)abab

t21%pos0 = posst
call grc0(4,4,3,3,1,1,1,t21%pos0,posst,t21%d,t21%i)
t22%pos0 = posst
call grc0(4,4,4,4,2,2,1,t22%pos0,posst,t22%d,t22%i)
t23%pos0 = posst
call grc0(4,0,3,4,1,2,1,t23%pos0,posst,t23%d,t23%i)

!1.7 maps for W0
!    W01 - <mn||ij>aaaa
!    W02 - <mn||ij>bbbb
!    W03 - <mn||ij>abab

w01%pos0 = posst
call grc0(4,4,1,1,1,1,1,w01%pos0,posst,w01%d,w01%i)
w02%pos0 = posst
call grc0(4,4,2,2,2,2,1,w02%pos0,posst,w02%d,w02%i)
w03%pos0 = posst
call grc0(4,0,1,2,1,2,1,w03%pos0,posst,w03%d,w03%i)

!1.8 maps for W1
!    W11 - <ie||mn>aaaa
!    W12 - <ie||mn>bbbb
!    W13 - <ie||mn>abab
!    W14 - <ie||mn>baab

w11%pos0 = posst
call grc0(4,3,1,3,1,1,1,w11%pos0,posst,w11%d,w11%i)
w12%pos0 = posst
call grc0(4,3,2,4,2,2,1,w12%pos0,posst,w12%d,w12%i)
w13%pos0 = posst
call grc0(4,0,1,4,1,2,1,w13%pos0,posst,w13%d,w13%i)
w14%pos0 = posst
call grc0(4,0,2,3,1,2,1,w14%pos0,posst,w14%d,w14%i)

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

v1%pos0 = posst
posst = posst+lengthv
v2%pos0 = posst
posst = posst+lengthv
v3%pos0 = posst
posst = posst+lengthv
v4%pos0 = posst
posst = posst+lengthv
lenv = lengthv

!2.2 M - files

m1%pos0 = posst
posst = posst+lengthm
m2%pos0 = posst
posst = posst+lengthm
m3%pos0 = posst
posst = posst+lengthm
m4%pos0 = posst
posst = posst+lengthm

!2.3 H - files

h1%pos0 = posst
posst = posst+lengthh
h2%pos0 = posst
posst = posst+lengthh
h3%pos0 = posst
posst = posst+lengthh
h4%pos0 = posst
posst = posst+lengthh

!2.4 N,P - files

n%pos0 = posst
posst = posst+lengthn
p%pos0 = posst
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
